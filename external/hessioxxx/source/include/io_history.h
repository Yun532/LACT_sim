/* ============================================================================

   Copyright (C) 1997, 2001, 2007, 2010, 2021  Konrad Bernloehr

   This file is part of the eventio/hessio library.

   The eventio/hessio library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this library. If not, see <http://www.gnu.org/licenses/>.

============================================================================ */

/** @file io_history.h
 *  @short Record history of configuration settings/commands.
 *
 *  @author  Konrad Bernloehr 
 *  @date 1997 to 2021
 */

#ifndef IO_HISTORY_H__LOADED
#define IO_HISTORY_H__LOADED 1

#ifndef INITIAL_H__LOADED
#include "initial.h"
#endif

#define IO_TYPE_HISTORY 70
#define IO_TYPE_CMD_HIST 71
#define IO_TYPE_CFG_HIST 72
#define IO_TYPE_METAPARAM 75

#ifdef __cplusplus
extern "C" {
#endif

/** Use to build a linked list of configuration history. */
struct history_struct
{
   char *text;  /** Configuration test */
   time_t time; /** Time when the configuration was entered */
   struct history_struct *next; /* Next element */
};
typedef struct history_struct HSTRUCT;

struct history_container_struct
{
   int id;                /**< Has always been 1 */
   HSTRUCT *cmdline;      /**< Prior commands executed for current data */
   HSTRUCT *cfg_global;   /**< Global (or rather: not telescope-specific) configuration */
   HSTRUCT **cfg_tel;     /**< One linked list per telescope, if recognized */
   size_t ntel;           /**< The number of telescopes for which cfg_tel was accumulated */
};
typedef struct history_container_struct HistoryContainer;

/* io_history.c */
int clear_hstruct(HSTRUCT *h);
int clear_histcont(HistoryContainer *hc);
int push_command_history (int argc, char **argv);
int push_config_history (const char *line, int noreplace);
int write_history (long id, IO_BUFFER *iobuf);
int read_history(IO_BUFFER *iobuf, HistoryContainer *hc);
int write_config_history (const char *htext, long htime, 
    long id, IO_BUFFER *iobuf);
int list_history (IO_BUFFER *iobuf, FILE *file);
int print_history (IO_BUFFER *iobuf);

/** A history meta parameter item consists of a parameter name and its
   text-formatted value, both as text strings. */
struct meta_param_item
{
   char *name;   /**< Parameter name */
   char *value;  /**< Parameter value in text representation */
   struct meta_param_item *next; /**< Pointer to next element in linked list or NULL for the last one. */
};
typedef struct meta_param_item MetaParamItem;

/** The linked MetaParamItem list for one ID (-1=global, detector ID otherwise)
   are registered under one list starting point. */
struct meta_param_list
{
   long ident;            /**< ID number (-1 for global parameter, 0 usually unused, >0 for detector ID) */
   MetaParamItem *first;  /**< Pointer to start of the linked list or NULL. */
};
typedef struct meta_param_list MetaParamList;

void set_metaparam_id (MetaParamList *lst, long id);
int fill_metaparam (MetaParamList *lst, const char **names, const char **values, size_t npar, long id);
int add_metaparam (MetaParamList *lst, const char *name, const char *value);
int clear_metaparam (MetaParamList *lst);
int write_metaparam (IO_BUFFER *iobuf, const MetaParamList *lst);
int read_metaparam (IO_BUFFER *iobuf, MetaParamList *lst);
int print_metaparam (IO_BUFFER *iobuf);
int show_metaparam (const MetaParamList *lst);
const char *search_metaparam (const MetaParamList *lst, const char *name);

#define WITH_METAPARAM_HISTORY 1

#ifdef __cplusplus
}
#endif

#endif
