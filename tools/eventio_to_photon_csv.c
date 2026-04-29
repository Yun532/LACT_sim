/*
 * Convert CORSIKA/EventIO photon bunch files to LACT_sim PhotonCsv.
 *
 * This is a thin adapter around the hessioxxx/eventio library. It does not
 * reimplement the EventIO binary format; hessioxxx owns block parsing,
 * compressed-file handling, and compact/long/3D photon-bunch decoding.
 */

#include "initial.h"
#include "io_basic.h"
#include "mc_tel.h"
#include "fileopen.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct ConverterOptions {
    const char *input_path;
    const char *output_path;
    const char *event_id_mode;
    int filter_event_id;
    int filter_array_id;
    int filter_telescope_id;
    int has_filter_event_id;
    int has_filter_array_id;
    int has_filter_telescope_id;
} ConverterOptions;

typedef struct ConverterStats {
    long blocks_seen;
    long photon_blocks_seen;
    long photon3d_blocks_seen;
    long rows_written;
    double multiplicity_sum;
} ConverterStats;

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s input.eventio.zst output.csv [options]\n"
        "\n"
        "Options:\n"
        "  --event-id-mode event|event_array100   Default: event\n"
        "  --filter-event-id N                    Keep only one output event id\n"
        "  --filter-array-id N                    Keep only one CORSIKA array block\n"
        "  --filter-telescope-id N                Keep only one telescope id\n"
        "\n"
        "Notes:\n"
        "  - hessioxxx returns photon positions in cm; PhotonCsv writes meters.\n"
        "  - event_array100 writes event_id = shower_event * 100 + array_id,\n"
        "    matching the ROOT runid convention used in the current validation file.\n",
        prog);
}

static int parse_int_arg(const char *text, int *out) {
    char *end = NULL;
    long value = strtol(text, &end, 10);
    if (end == text || *end != '\0') {
        return -1;
    }
    *out = (int)value;
    return 0;
}

static int parse_args(int argc, char **argv, ConverterOptions *opt) {
    memset(opt, 0, sizeof(*opt));
    opt->event_id_mode = "event";

    if (argc < 3) {
        usage(argv[0]);
        return -1;
    }
    opt->input_path = argv[1];
    opt->output_path = argv[2];

    for (int i = 3; i < argc; ++i) {
        if (strcmp(argv[i], "--event-id-mode") == 0 && i + 1 < argc) {
            opt->event_id_mode = argv[++i];
            if (strcmp(opt->event_id_mode, "event") != 0 &&
                strcmp(opt->event_id_mode, "event_array100") != 0) {
                fprintf(stderr, "unknown --event-id-mode: %s\n", opt->event_id_mode);
                return -1;
            }
        } else if (strcmp(argv[i], "--filter-event-id") == 0 && i + 1 < argc) {
            if (parse_int_arg(argv[++i], &opt->filter_event_id) != 0) {
                fprintf(stderr, "invalid --filter-event-id\n");
                return -1;
            }
            opt->has_filter_event_id = 1;
        } else if (strcmp(argv[i], "--filter-array-id") == 0 && i + 1 < argc) {
            if (parse_int_arg(argv[++i], &opt->filter_array_id) != 0) {
                fprintf(stderr, "invalid --filter-array-id\n");
                return -1;
            }
            opt->has_filter_array_id = 1;
        } else if (strcmp(argv[i], "--filter-telescope-id") == 0 && i + 1 < argc) {
            if (parse_int_arg(argv[++i], &opt->filter_telescope_id) != 0) {
                fprintf(stderr, "invalid --filter-telescope-id\n");
                return -1;
            }
            opt->has_filter_telescope_id = 1;
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            usage(argv[0]);
            return -1;
        } else {
            fprintf(stderr, "unknown option or missing value: %s\n", argv[i]);
            usage(argv[0]);
            return -1;
        }
    }
    return 0;
}

static int output_event_id(int shower_event_id, int array_id, const ConverterOptions *opt) {
    if (strcmp(opt->event_id_mode, "event_array100") == 0) {
        return shower_event_id * 100 + array_id;
    }
    return shower_event_id;
}

static int should_keep(int event_id, int array_id, int telescope_id, const ConverterOptions *opt) {
    if (opt->has_filter_event_id && event_id != opt->filter_event_id) {
        return 0;
    }
    if (opt->has_filter_array_id && array_id != opt->filter_array_id) {
        return 0;
    }
    if (opt->has_filter_telescope_id && telescope_id != opt->filter_telescope_id) {
        return 0;
    }
    return 1;
}

static double downward_dir_z(double cx, double cy) {
    double z2 = 1.0 - cx * cx - cy * cy;
    if (z2 < 0.0) {
        z2 = 0.0;
    }
    return -sqrt(z2);
}

static int read_current_event_id(IO_BUFFER *iobuf, int item_type, int *current_event_id) {
    real data[273];
    if (read_tel_block(iobuf, item_type, data, 273) < 0) {
        return -1;
    }
    if (item_type == IO_TYPE_MC_EVTH) {
        *current_event_id = (int)Nint(data[1]);
    }
    return 0;
}

static int write_tel_photon_block(IO_BUFFER *iobuf, FILE *out, const ConverterOptions *opt,
                                  int current_event_id, ConverterStats *stats) {
    int array_id = 0;
    int telescope_id = 0;
    int nbunches = 0;
    double photons = 0.0;
    int rc = read_tel_photons(iobuf, 0, &array_id, &telescope_id, &photons, NULL, &nbunches);
    if (rc != -10) {
        return rc;
    }

    stats->photon_blocks_seen++;
    int event_id = output_event_id(current_event_id, array_id, opt);
    if (!should_keep(event_id, array_id, telescope_id, opt)) {
        struct bunch *skip = (struct bunch *)malloc((size_t)nbunches * sizeof(struct bunch));
        if (skip == NULL && nbunches > 0) {
            fprintf(stderr, "failed to allocate skip buffer for %d bunches\n", nbunches);
            return -1;
        }
        rc = read_tel_photons(iobuf, nbunches, &array_id, &telescope_id, &photons, skip, &nbunches);
        free(skip);
        return rc;
    }

    struct bunch *bunches = (struct bunch *)malloc((size_t)nbunches * sizeof(struct bunch));
    if (bunches == NULL && nbunches > 0) {
        fprintf(stderr, "failed to allocate buffer for %d bunches\n", nbunches);
        return -1;
    }
    rc = read_tel_photons(iobuf, nbunches, &array_id, &telescope_id, &photons, bunches, &nbunches);
    if (rc < 0) {
        free(bunches);
        return rc;
    }

    event_id = output_event_id(current_event_id, array_id, opt);
    for (int i = 0; i < nbunches; ++i) {
        const struct bunch *b = &bunches[i];
        fprintf(out, "%.9g,%.9g,0,%.9g,%.9g,%.12g,%.9g,%.9g,1,%.9g,%d,%d\n",
                b->x * 0.01, b->y * 0.01,
                b->cx, b->cy, downward_dir_z(b->cx, b->cy),
                b->ctime, b->lambda, b->photons, event_id, telescope_id);
        stats->rows_written++;
        stats->multiplicity_sum += b->photons;
    }
    free(bunches);
    return 0;
}

static int write_tel_photon3d_block(IO_BUFFER *iobuf, FILE *out, const ConverterOptions *opt,
                                    int current_event_id, ConverterStats *stats) {
    int array_id = 0;
    int telescope_id = 0;
    int nbunches = 0;
    double photons = 0.0;
    int rc = read_tel_photons3d(iobuf, 0, &array_id, &telescope_id, &photons, NULL, &nbunches);
    if (rc != -10) {
        return rc;
    }

    stats->photon3d_blocks_seen++;
    int event_id = output_event_id(current_event_id, array_id, opt);
    if (!should_keep(event_id, array_id, telescope_id, opt)) {
        struct bunch3d *skip = (struct bunch3d *)malloc((size_t)nbunches * sizeof(struct bunch3d));
        if (skip == NULL && nbunches > 0) {
            fprintf(stderr, "failed to allocate 3D skip buffer for %d bunches\n", nbunches);
            return -1;
        }
        rc = read_tel_photons3d(iobuf, nbunches, &array_id, &telescope_id, &photons, skip, &nbunches);
        free(skip);
        return rc;
    }

    struct bunch3d *bunches = (struct bunch3d *)malloc((size_t)nbunches * sizeof(struct bunch3d));
    if (bunches == NULL && nbunches > 0) {
        fprintf(stderr, "failed to allocate 3D buffer for %d bunches\n", nbunches);
        return -1;
    }
    rc = read_tel_photons3d(iobuf, nbunches, &array_id, &telescope_id, &photons, bunches, &nbunches);
    if (rc < 0) {
        free(bunches);
        return rc;
    }

    event_id = output_event_id(current_event_id, array_id, opt);
    for (int i = 0; i < nbunches; ++i) {
        const struct bunch3d *b = &bunches[i];
        fprintf(out, "%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,1,%.9g,%d,%d\n",
                b->x * 0.01, b->y * 0.01, b->z * 0.01,
                b->cx, b->cy, b->cz, b->ctime, b->lambda,
                b->photons, event_id, telescope_id);
        stats->rows_written++;
        stats->multiplicity_sum += b->photons;
    }
    free(bunches);
    return 0;
}

static int handle_tel_array(IO_BUFFER *iobuf, FILE *out, const ConverterOptions *opt,
                            int current_event_id, ConverterStats *stats) {
    IO_ITEM_HEADER array_header;
    int array_id = 0;
    int type = 0;
    int rc = begin_read_tel_array(iobuf, &array_header, &array_id);
    if (rc < 0) {
        return rc;
    }
    while ((type = next_subitem_type(iobuf)) > 0) {
        if (type == IO_TYPE_MC_PHOTONS) {
            rc = write_tel_photon_block(iobuf, out, opt, current_event_id, stats);
        } else if (type == IO_TYPE_MC_PHOTONS3D) {
            rc = write_tel_photon3d_block(iobuf, out, opt, current_event_id, stats);
        } else {
            rc = skip_subitem(iobuf);
        }
        if (rc < 0) {
            get_item_end(iobuf, &array_header);
            return rc;
        }
    }
    return end_read_tel_array(iobuf, &array_header);
}

int main(int argc, char **argv) {
    ConverterOptions opt;
    ConverterStats stats;
    IO_BUFFER *iobuf = NULL;
    IO_ITEM_HEADER item_header;
    FILE *out = NULL;
    int current_event_id = 0;
    int rc = 0;

    if (parse_args(argc, argv, &opt) != 0) {
        return 2;
    }
    memset(&stats, 0, sizeof(stats));

    iobuf = allocate_io_buffer(5000000L);
    if (iobuf == NULL) {
        fprintf(stderr, "cannot allocate EventIO buffer\n");
        return 1;
    }
    if (iobuf->max_length < 1000000000L) {
        iobuf->max_length = 1000000000L;
    }

    iobuf->input_file = fileopen(opt.input_path, READ_BINARY);
    if (iobuf->input_file == NULL) {
        fprintf(stderr, "cannot open input %s: %s\n", opt.input_path, strerror(errno));
        free_io_buffer(iobuf);
        return 1;
    }

    out = fopen(opt.output_path, "w");
    if (out == NULL) {
        fprintf(stderr, "cannot open output %s: %s\n", opt.output_path, strerror(errno));
        fileclose(iobuf->input_file);
        free_io_buffer(iobuf);
        return 1;
    }
    fprintf(out, "x_m,y_m,z_m,dir_x,dir_y,dir_z,time_ns,wavelength_nm,weight,multiplicity,event_id,telescope_id\n");

    while (find_io_block(iobuf, &item_header) == 0) {
        if (read_io_block(iobuf, &item_header) != 0) {
            break;
        }
        stats.blocks_seen++;
        switch ((int)item_header.type) {
            case IO_TYPE_MC_EVTH:
            case IO_TYPE_MC_RUNH:
            case IO_TYPE_MC_EVTE:
            case IO_TYPE_MC_RUNE:
                rc = read_current_event_id(iobuf, (int)item_header.type, &current_event_id);
                break;
            case IO_TYPE_MC_TELARRAY:
                rc = handle_tel_array(iobuf, out, &opt, current_event_id, &stats);
                break;
            case IO_TYPE_MC_PHOTONS:
                rc = write_tel_photon_block(iobuf, out, &opt, current_event_id, &stats);
                break;
            case IO_TYPE_MC_PHOTONS3D:
                rc = write_tel_photon3d_block(iobuf, out, &opt, current_event_id, &stats);
                break;
            default:
                rc = 0;
                break;
        }
        if (rc < 0) {
            fprintf(stderr, "error while reading EventIO block type %lu (rc=%d)\n",
                    item_header.type, rc);
            fclose(out);
            fileclose(iobuf->input_file);
            free_io_buffer(iobuf);
            return 1;
        }
    }

    fclose(out);
    fileclose(iobuf->input_file);
    free_io_buffer(iobuf);

    fprintf(stderr,
            "eventio_to_photon_csv summary:\n"
            "  input              : %s\n"
            "  output             : %s\n"
            "  blocks_seen        : %ld\n"
            "  photon_blocks      : %ld\n"
            "  photon3d_blocks    : %ld\n"
            "  rows_written       : %ld\n"
            "  multiplicity_sum   : %.12g\n",
            opt.input_path, opt.output_path, stats.blocks_seen,
            stats.photon_blocks_seen, stats.photon3d_blocks_seen,
            stats.rows_written, stats.multiplicity_sum);

    return 0;
}
