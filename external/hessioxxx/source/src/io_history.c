/* ============================================================================

   Copyright (C) 1997, 2001, 2005, 2007, 2009, 2010, 2014, 2021, 2022  Konrad Bernloehr

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

/** @file io_history.c
 *  @short Record history of configuration settings/commands.
 *
 *  This code has not been adapted for multi-threading.
 *
 *  @author  Konrad Bernloehr 
 *  @date 1997 to 2022
 */
/* ================================================================ */

#include "initial.h"
#include "io_basic.h"
#include "io_history.h"
// #include "history.h"
#include "current.h"
#include "hconfig.h"

static char *cmdline = NULL; /**< A copy of the program's command line. */
static time_t cmdtime;       /**< The time when the program was started. */
static HSTRUCT *configs = NULL; /**< Start of configuration history. */

static void listtime (time_t t,FILE * f);

/* ------------------------ push_command_history --------------------- */
/**
 *    @short Save the command line for later output in a history block.
 *
 *    @param  argc  Number of command line arguments (incl. command)
 *    @param  argv  Pointers to argument text strings
 *
 *    @return  0 (o.k.), -1 (invalid argument or no memory)
 *
 */

int push_command_history (int argc, char **argv)
{
   int i;
   size_t len;

   for (i=0, len=0; i<argc && argv[i] != (char *) NULL; i++)
      len += strlen(argv[i])+1;
   if ( len < 1 )
      return -1;

   if ( cmdline != (char *) NULL )
      free((void *)cmdline);
   if ( (cmdline = (char *) malloc(len)) == (char *) NULL )
      return -1;
   strcpy(cmdline,argv[0]);
   for (i=1; i<argc && argv[i] != (char *) NULL; i++)
   {
      strcat(cmdline," ");
      strcat(cmdline,argv[i]);
   }
   cmdtime = current_time();

   return 0;
}

/* -------------------------- push_config_history ---------------------- */
/**
 *    @short Save a line of configuration text for later output in a history block.
 *    If any configuration text for the same keyword was present before,
 *    it may be replaced by the new text, even if the new setting does not
 *    replace all old settings.
 *
 *    @param  line    Configuration text line.
 *    @param  replace Replace old text for same keyword (1) or not (0).
 *
 *    @return 0 (o.k.),  -1 (memory allocation failed)
 *
 */

int push_config_history (const char *line, int noreplace)
{
   char word[32], word2[32];
   char *text;
   int ipos;
   int replace;
   HSTRUCT *hp, *htmp;
   if ( noreplace )
      replace = 0;
   else
      replace = 1;

   if ( line == (char *) NULL )
      return -1;

   ipos = 0;
   if ( getword(line,&ipos,word,sizeof(word)-1,' ','%') <= 0 )
      return 0;

   if ( (text = (char*) malloc(strlen(line)+1)) == (char *) NULL )
      return -1;
   strcpy(text,line);

   if ( configs == (HSTRUCT *) NULL )
   {
      
      if ( (configs = (HSTRUCT *) calloc(1,sizeof(HSTRUCT))) == 
           (HSTRUCT *) NULL )
      {
         free((void *)text);
         return -1;
      }
      configs->text = text;
      configs->time = current_time();
      configs->next = (HSTRUCT *) NULL;
   }
   else
   {
      for (hp=configs; hp->next != (HSTRUCT *) NULL; hp=hp->next)
      {
         if ( replace )
         {
            ipos = 0;
            if ( getword(hp->text,&ipos,word2,sizeof(word2)-1,' ','%') <= 0 )
               continue;
            /* If a line for the same keyword exists, remove it. */
            if ( strcmp(word,word2) == 0 )
            {
               free((void *)hp->text);
               hp->text = text;
               hp->time = current_time();
               return 0;
            }
         }
      }
      if ( replace )
      {
         ipos = 0;
         if ( getword(hp->text,&ipos,word2,sizeof(word2)-1,' ','%') > 0 )
            /* If a line for the same keyword exists, remove it. */
            if ( strcmp(word,word2) == 0 )
            {
               free((void *)hp->text);
               hp->text = text;
               hp->time = current_time();
               return 0;
            }
      }

      if ( (htmp = (HSTRUCT *) calloc(1,sizeof(HSTRUCT))) == 
           (HSTRUCT *) NULL )
      {
         free((void *)text);
         return -1;
      }
      hp->next = htmp;
      htmp->text = text;
      htmp->time = current_time();
      htmp->next = (HSTRUCT *) NULL;
   }

   return 0;
}

/* -------------------------- clear_hstruct ---------------------- */
/**
 *  @short Clear and free all elements along one linked list of history elements.
 *
 *  @param  h Start of linked list of HSTRUCTs.
 *  @return Number of elements released.
 */

int clear_hstruct(HSTRUCT *h)
{
   HSTRUCT *next = h;
   int n = 0;
   while ( next != NULL )
   {
      HSTRUCT *h1 = next;
      next = h1->next;
      if ( h1->text != NULL )
      {
         free(h1->text);
         h1->text = NULL;
      }
      h1->time = 0;
      h1->next = NULL;
      free(h1);
      n++;
   }
   return n;
}

/* -------------------------- clear_histcont ---------------------- */
/**
 *  @short Clear and free all linked list of history elements in a history container.
 *  The container itself will only contain null pointers afterwards but no
 *  attempts are made to free() it as well - it could just as well come from the stack.
 *
 *  @param  hc Pointer to history container.
 *  @return Number of elements released.
 */

int clear_histcont(HistoryContainer *hc)
{
   int n = 0;
   size_t itel = 0;

   if ( hc == NULL )
      return 0;

   if ( hc->cmdline != NULL )
   {
      n += clear_hstruct(hc->cmdline);
      hc->cmdline = NULL; /* All elements already free'ed */
   }
   if ( hc->cfg_global != NULL )
   {
      n += clear_hstruct(hc->cfg_global);
      hc->cfg_global = NULL; /* All elements already free'ed */
   }
   if ( hc->cfg_tel != NULL )
   {
      for ( itel=0; itel<hc->ntel; itel++ )
      {
         if ( hc->cfg_tel[itel] != NULL )
         {
            n += clear_hstruct(hc->cfg_tel[itel]);
            hc->cfg_tel[itel] = NULL; /* All elements already free'ed */
         }
      }
      free(hc->cfg_tel);
      hc->cfg_tel = NULL;
   }
   hc->id = 0;
   /* We leave the empty container, not making any assumptions on it
      coming from the stack or heap. */
   
   return n;
}

/* -------------------------- write_history ---------------------- */
/**
 *    @short Write the block of accumulated history lines (command line and
 *    configuration lines) to an I/O buffer.
 *
 *    @param   id    Identifier (detector number)
 *    @param   iobuf I/O buffer descriptor
 *
 *    @return  0 (o.k.), <0 (error with I/O buffer)
 *
 */

int write_history (long id, IO_BUFFER *iobuf)
{
   IO_ITEM_HEADER item_header, sub_item_header;
   HSTRUCT *hp;

   item_header.type = IO_TYPE_HISTORY;  /* History block is type 70 */
   item_header.version = 1;             /* Version 1 */
   item_header.ident = id;

   if ( put_item_begin(iobuf,&item_header) < 0 )
      return -1;

   if ( cmdline != (char *) NULL )
   {
      sub_item_header.type = IO_TYPE_CMD_HIST;  /* Command line history is type 71 */
      sub_item_header.version = 1;              /* Version 1 */
      sub_item_header.ident = id;
      if ( put_item_begin(iobuf,&sub_item_header) < 0 )
         return -1;
      put_long(cmdtime,iobuf);
      put_string(cmdline,iobuf);
      put_item_end(iobuf,&sub_item_header);
   }

   if ( configs != (HSTRUCT *) NULL )
      for (hp=configs; hp!=(HSTRUCT *) NULL; hp=hp->next)
         if ( write_config_history(hp->text,hp->time,id,iobuf) < 0 )
            break;

   return(put_item_end(iobuf,&item_header));
}

/* -------------------------- read_history ---------------------- */
/**
 *    @short Read the block of accumulated history lines (command line and
 *    configuration lines) from an I/O buffer and try to split the
 *    configuration history by telescope, if the transition can be found.
 *
 *    Note that in extracted/split/merged sim_telarray data the correspondence to
 *    current telescope IDs is lost. For other programs, no
 *    telescope transition is recognized, i.e. all configuration is global.
 *
 *    @param   iobuf I/O buffer descriptor
 *    @param   hc    Pointer to history container.
 *
 *    @return  0 (o.k.), <0 (error with I/O buffer)
 *
 */

int read_history (IO_BUFFER *iobuf, HistoryContainer *hc)
{
   IO_ITEM_HEADER item_header, sub_item_header;
   int type = 0, rc = 0, itel = 0;
   char line[10240];
   time_t t;
   HSTRUCT *pos = NULL;
   
   if ( iobuf == (IO_BUFFER *) NULL || hc == NULL )
      return -1;

   /* Make sure the container is clean before we fill anything in. */
   clear_histcont(hc);

   item_header.type = IO_TYPE_HISTORY;  /* Data type */
   if ( (rc = get_item_begin(iobuf,&item_header)) < 0 )
      return rc;
   if ( item_header.version > 1 )
   {
      fprintf(stderr,"Unsupported history version: %d.\n",
         item_header.version);
      get_item_end(iobuf,&item_header);
      return -1;
   }
   
   hc->id = item_header.ident; /* Has always been '1' so far. */
   
   while ( (type = next_subitem_type(iobuf)) > 0 )
   {
      if ( type == IO_TYPE_CMD_HIST )
      {
         sub_item_header.type = IO_TYPE_CMD_HIST;  /* Command line history is type 71 */
         if ( (rc = get_item_begin(iobuf,&sub_item_header)) < 0 )
         {
            get_item_end(iobuf,&item_header);
            return rc;
         }
         if ( sub_item_header.version > 1 )
         {
            fprintf(stderr,"Unsupported command-line history version: %d.\n",
               sub_item_header.version);
            get_item_end(iobuf,&sub_item_header);
            get_item_end(iobuf,&item_header);
            return -1;
         }
         t = get_long(iobuf);
         get_string(line,sizeof(line)-1,iobuf);
         if ( (hc->cmdline = (HSTRUCT *) calloc(1,sizeof(HSTRUCT))) != NULL )
         {
            hc->cmdline->time = t;
            hc->cmdline->text = strdup(line);
         }
      }
      else if ( type == IO_TYPE_CFG_HIST )
      {
         sub_item_header.type = IO_TYPE_CFG_HIST;  /* Config history is type 72 */
         if ( (rc = get_item_begin(iobuf,&sub_item_header)) < 0 )
         {
            get_item_end(iobuf,&item_header);
            return rc;
         }
         if ( sub_item_header.version > 1 )
         {
            fprintf(stderr,"Unsupported config history version: %d.\n",
               sub_item_header.version);
            get_item_end(iobuf,&sub_item_header);
            get_item_end(iobuf,&item_header);
            return -1;
         }
         t = get_long(iobuf);
         get_string(line,sizeof(line)-1,iobuf);
         if ( pos == NULL ) /* First one is always starting global config */
         {
            if ( (hc->cfg_global = (HSTRUCT *) malloc(sizeof(HSTRUCT))) == NULL )
               return -4;

            pos = hc->cfg_global;
            pos->time = t;
            pos->text = strdup(line);
            pos->next = NULL;
         }
         else if ( strncmp(line,"# Telescope-specific configuration follows",40) == 0 )
         {
            /* Moving to the first or next telescope */
            if ( hc->ntel <= 0 )
            {
               if ( (hc->cfg_tel = (HSTRUCT **) malloc(sizeof(HSTRUCT *))) == NULL )
                  return -4;
               hc->ntel = 1;
               itel = 0;
            }
            else
            {
               HSTRUCT **o = hc->cfg_tel;
               if ( (hc->cfg_tel = (HSTRUCT **) realloc(hc->cfg_tel,(hc->ntel+1)*sizeof(HSTRUCT *))) == NULL )
               {
                  hc->cfg_tel = o;
                  return -4;
               }
               itel = hc->ntel;
               hc->ntel++;
            }
            if ( (hc->cfg_tel[itel] = (HSTRUCT *) malloc(sizeof(HSTRUCT))) == NULL )
               return -4;
            pos = hc->cfg_tel[itel];
            pos->time = t;
            pos->text = strdup(line);
            pos->next = NULL;
         }
         else
         {
            /* Another configuration item in the current section */
            if ( (pos->next = (HSTRUCT *) malloc(sizeof(HSTRUCT))) == NULL )
               return -4;
            pos = pos->next;
            pos->time = t;
            pos->text = strdup(line);
            pos->next = NULL;
         }
      }
      else
      {
         fprintf(stderr,"Unsupported history sub-block type: %d.\n", type);
         get_item_end(iobuf,&item_header);
         return -1;
      }
   }
   
   return (get_item_end(iobuf,&item_header));
}

/* ------------------------ write_config_history -------------------- */
/**
 *    @short Write a configuration history line to an I/O buffer.
 *
 *    @param   htext  Text of configuration line
 *    @param   htime  Time when the configuration was set.
 *    @param   id     Identifier (detector number)
 *    @param   iobuf  I/O buffer descriptor
 *
 *    @return  0 (o.k.), <0 (error with I/O buffer)
 *
 */

int write_config_history (const char *htext, long htime, long id, IO_BUFFER *iobuf)
{
   IO_ITEM_HEADER item_header;

   item_header.type = IO_TYPE_CFG_HIST;  /* Configuration line is type 72 */
   item_header.version = 1;              /* Version 1 */
   item_header.ident = id;

   if ( put_item_begin(iobuf,&item_header) < 0 )
      return -1;
   put_long(htime,iobuf);
   put_string(htext,iobuf);

   return(put_item_end(iobuf,&item_header));
}

/* ------------------------- listtime ----------------------- */

static void listtime (time_t t, FILE *f)
{
   struct tm *tmp;
   struct tm tmbuf;
   static int once = 1;

   if ( once )
   {
      tzset(); /* Not needed for localtime() but for localtime_r() */
      once = 0;
   }

   if ( (tmp = localtime_r(&t,&tmbuf)) == (struct tm *) NULL )
      return;

#if defined _POSIX_VERSION && _POSIX_VERSION >= 200112L
   /* POSIX.1 (or BSD) has extern char *tzname[2] for std/dst names of timezone */
   const char *tzn = (tmp->tm_isdst ? tzname[1] : tzname[0]);
#else
   const char *tzn = "UTC?";
#endif
/* There may be other fields in struct tm holding the actual timezone name */
#ifdef __GLIBC__
# ifdef	__USE_MISC
   tzn = tmp->tm_zone;
# else
   tzn = tmp->__tm_zone;
# endif
#endif

   fprintf(f,"%04d-%02d-%02d %02d:%02d:%02d %s",
      tmp->tm_year+1900, tmp->tm_mon+1, tmp->tm_mday, 
      tmp->tm_hour, tmp->tm_min, tmp->tm_sec, tzn);
}

/* ------------------------ list_history ----_------------------------- */
/**
 *    @short List history block contents on standard output or other file.
 *
 *    @param  iobuf  I/O buffer descriptor
 *    @param  file   Optional open output stream (NULL -> stdout).
 *
 *    @return 0 (o.k.), <0 (error with I/O buffer)
 *
 */

int list_history (IO_BUFFER *iobuf, FILE *file)
{
   IO_ITEM_HEADER item_header, sub_item_header;
   int rc;
   char line[10240];
   char *s;
   time_t t;
   
   if ( file == NULL )
      file = stdout;
  
   if ( iobuf == (IO_BUFFER *) NULL )
      return -1;
      
   item_header.type = IO_TYPE_HISTORY;
   if ( (rc = get_item_begin(iobuf,&item_header)) < 0 )
      return rc;
   if ( item_header.version != 1 )
   {
      Warning("Wrong version number of history item to be read.");
      return -1;
   }
   if ( item_header.ident > 0 )
      fprintf(file,"\nHistory block (ident %ld):\n",item_header.ident);
   else
      fprintf(file,"\nHistory block:\n");
   sub_item_header.type = IO_TYPE_CMD_HIST;
   if ( search_sub_item(iobuf,&item_header,&sub_item_header) == 0 )
   {
      if ( (rc = get_item_begin(iobuf,&sub_item_header)) < 0 )
      {
         get_item_end(iobuf,&item_header);
         return rc;
      }
      if ( sub_item_header.version != 1 )
      {
         Warning("Wrong version number of command line history item to be read.");
         get_item_end(iobuf,&sub_item_header);
         get_item_end(iobuf,&item_header);
         return -1;
      }
      t = get_long(iobuf);
      get_string(line,sizeof(line)-1,iobuf);
      fprintf(file,"   Command line (dated ");
      listtime(t,file);
      if ( sub_item_header.ident > 0 )
         fprintf(file,", ident %ld):\n      %s\n",sub_item_header.ident,line);
      else
         fprintf(file,"):\n      %s\n",line);
      get_item_end(iobuf,&sub_item_header);
   }
   rewind_item(iobuf,&item_header);
   fprintf(file,"   Configuration data:\n");
   for (;;)
   {
      sub_item_header.type = IO_TYPE_CFG_HIST;
      if ( search_sub_item(iobuf,&item_header,&sub_item_header) != 0 )
         break;
      if ( (rc = get_item_begin(iobuf,&sub_item_header)) < 0 )
      {
         get_item_end(iobuf,&item_header);
         return rc;
      }
      if ( sub_item_header.type != IO_TYPE_CFG_HIST )
      {
         Warning("Wrong item type instead of configuration history item");
         return -1;
      }
      if ( sub_item_header.version != 1 )
      {
         Warning("Wrong version number of configuration history item to be read.");
         get_item_end(iobuf,&sub_item_header);
         get_item_end(iobuf,&item_header);
         return -1;
      }
      *line = '\0';
      t = get_long(iobuf);
      get_string(line,sizeof(line)-1,iobuf);
      fprintf(file,"      ");
      for (s=line; *s!='\0'; s++)
      {
         fputc(*s,file);
         if ( *s=='\n' && *(s+1)!='\0' )
            fputs("      ",file);
      }
      if ( s!=line && *(s-1) != '\n' )
         fputc('\n',file);

      get_item_end(iobuf,&sub_item_header);
   }
   fputs("End of history block\n",file);

   return (get_item_end(iobuf,&item_header));
}

/* ------------------------ print_history ----_------------------------- */
/**
 *    @short Stub for list_history(), to stdout only.
 *
 *    @param  iobuf  I/O buffer descriptor
 *
 *    @return 0 (o.k.), <0 (error with I/O buffer)
 *
 */

int print_history (IO_BUFFER *iobuf)
{
   return list_history(iobuf,stdout);
}

/* ========== Meta-parameter part of the program history data ================= */
/*    Preserve selected parameters and their values no matter if they have
 *    any direct purpose for running the program or only as a reminder
 *    and for descriptive purposes.
 */

void set_metaparam_id (MetaParamList *lst, long id) { if ( lst != NULL ) lst->ident = id; }

/* ------------------------- fill_metaparam --------------------- */
/**
 *  @short Fill in meta parameters to the linked list after (re-)initialising it.
 *
 *  @param lst The starting point of a linked list. Must exist but can be without items.
 *  @param names Array of names of the items to add. Array and all names must exist and be non-empty.
 *  @paran values Array of values of the items to add. Array and and values must exist (but can be empty).
 *
 *  @return 0 (OK), -1 or other non-zero (problem)
 */

int fill_metaparam (MetaParamList *lst, const char **names, const char **values, size_t npar, long id)
{
   MetaParamItem *item = NULL;
   size_t ipar;
   int rc = 0;

   if ( lst == (MetaParamList *) NULL || names == NULL || values == NULL )
      return -1;

   /* We always start with a fresh list */
   clear_metaparam(lst);
   lst->ident = id;

   if ( npar == 0 )
      return 0; /* No need to book the first item if nothing to be added */

   /* Check for any invalid entries: results in empty list */
   for (ipar=0; ipar<npar; ipar++)
   {
      if ( names[ipar] == NULL || names[ipar][0] == '\0' || values[ipar] == NULL )
      {
         rc = -1;
         break;
      }
   }

   if ( rc != 0 ) /* Any problems? */
      return rc;

   item = lst->first = (MetaParamItem *) calloc(1,sizeof(MetaParamItem));
   if ( item == NULL )
   {
      clear_metaparam(lst);
      return -1;
   }
   
   for (ipar=0; ipar+1<npar; ipar++)
   {
      item->name = strdup(names[ipar]); /* Note: not checking for duplicates */
      item->value = strdup(values[ipar]);
      item->next = (MetaParamItem *) calloc(1,sizeof(MetaParamItem));
      if ( item->name == NULL || item->value == NULL || item->next == NULL )
      {
         clear_metaparam(lst);
         return -1;
      }
      item = item->next;
   }
   /* The last item does not need to extend the list */
   ipar = npar-1;
   item->name = strdup(names[ipar]); /* Note: not checking for duplicates */
   item->value = strdup(values[ipar]);
   if ( item->name == NULL || item->value == NULL )
   {
      clear_metaparam(lst);
      return -1;
   }
   
   return 0;
}

/* -------------------------- add_metaparam --------------------- */
/**
 *  @short Add a meta parameter to the linked list or change a matching enty.
 *
 *  @param lst The starting point of a linked list. Must exist but can be without items.
 *  @param name Name of the item to add. Must exist and be non-empty.
 *  @param value Value of the item to add. Must exist.
 *
 *  @return 0 (OK), -1 or other non-zero (problem)
 */

int add_metaparam (MetaParamList *lst, const char *name, const char *value)
{
   MetaParamItem *item = NULL, *next = NULL;

   if ( lst == (MetaParamList *) NULL || name == NULL || value == NULL )
      return -1;
   if ( name[0] == '\0' ) /* Name must not be empty (but value can be) */
      return -1;

   /* No items filled in at all? */
   if ( lst->first == (MetaParamItem *) NULL )
   {
      /* Book the first item in the linked list and finish */
      if ( (lst->first = (MetaParamItem *) calloc(1,sizeof(MetaParamItem))) == NULL )
         return -1;
      item = lst->first;
      item->name = strdup(name);
      item->value = strdup(value);
      if ( item->name == NULL || item->value == NULL )
         return -1;
      return 0;
   }

   /* Go over all items, check for matches or find the end of it */
   for ( item = lst->first; item != NULL; item = next)
   {
      next = item->next;
      if ( item->name == NULL || item->value == NULL )
         continue;
      /* Does this item match the provided name? */
      if ( strcmp(item->name,name) == 0 )
      {
         /* Is it the same value again? */
         if ( strcmp(item->value,value) != 0 )
         {
            /* Only replacing if the value has changed */
            free(item->value);
            item->value = strdup(value);
         }
         return 0;
      }
      if ( next == NULL ) /* Break out of loop before we lose item pointer */
         break;
   }
   
   /* None of the available items matched the provided name */
   if ( item == (MetaParamItem *) NULL )
      return -1;

   /* Append a new item to the end of the linked list */
   if ( (item->next = (MetaParamItem *) calloc(1,sizeof(MetaParamItem))) == NULL )
      return -1;
   item = item->next;
   item->name = strdup(name);
   item->value = strdup(value);
   if ( item->name == NULL || item->value == NULL )
      return -1;
   
   return 0;
}

/* ----------------------- clear_metaparam ---------------- */
/** @short Clear a list of meta parameters. All pointers along
 *         the linked list must be dynamically allocated.
 *
 *  @param lst The starting point of a linked list. Must exist.
 *
 *  @return 0 (OK), -1 or other non-zero (problem)
 */

int clear_metaparam (MetaParamList *lst)
{
   MetaParamItem *item = NULL, *next = NULL;

   if ( lst == (MetaParamList *) NULL )
      return -1;
   
   for ( item = lst->first; item != NULL; item = next)
   {
      next = item->next;
      if ( item->name != NULL )
      {
         *(item->name) = '\0';
         free(item->name);
      }
      if ( item->value != NULL )
      {
         *(item->value) = '\0';
         free(item->value);
      }
      if ( next == NULL )
         break;
      item->next = NULL;
      free(item);
   }
   
   lst->ident = 0;
   
   return 0;
}

/* --------------------------- write_metaparam ------------------- */
/** @short Write a data block of meta parameters.
 *
 *  @param iobuf I/O buffer descriptor
 *  @param lst The starting point of a linked list. Must exist.
 *
 *  @return 0 (OK), -1 or other non-zero (problem)
 */

int write_metaparam (IO_BUFFER *iobuf, const MetaParamList *lst)
{
   MetaParamItem *item = NULL, *next = NULL;
   IO_ITEM_HEADER item_header;
   size_t n=0;

   item_header.type = IO_TYPE_METAPARAM;        /* Metadata block is type 75 */
#ifdef WITH_METAPARAM_CHECK
   item_header.version = 1;         /* V1 not intended for actual use */
#else
   item_header.version = 0;         /* Version 0 */
#endif
   item_header.ident = lst->ident;
   
   if ( lst == (const MetaParamList *) NULL || iobuf == (IO_BUFFER *) NULL )
      return -1;

   if ( put_item_begin(iobuf,&item_header) < 0 )
      return -1;

   /* Count the items that will be written */
   for ( item = lst->first; item != NULL; item = next)
   {
      next = item->next;
      if ( item->name == NULL )
         continue;
      n++;
   }
   
   /* Write down the number of items to be expected when reading */
   put_count(n, iobuf);

   /* Now we can write the items one by one */
   for ( item = lst->first; item != NULL; item = next)
   {
      next = item->next;
      if ( item->name == NULL )
         continue;
      put_var_string(item->name,iobuf);
      if ( item->value != NULL )
         put_var_string(item->value,iobuf);
      else
         put_var_string("",iobuf);
   }

   return(put_item_end(iobuf,&item_header));
}

/* --------------------------- read_metaparam ------------------- */
/** @short Read from a data block of meta parameters.
 *
 *  @param iobuf I/O buffer descriptor
 *  @param lst The starting point of a linked list. Must exist.
 *
 *  @return 0 (OK), -1 or other non-zero (problem)
 */

int read_metaparam (IO_BUFFER *iobuf, MetaParamList *lst)
{
   MetaParamItem *item = NULL;
   IO_ITEM_HEADER item_header;
   char item_name[1024], item_value[10240];
   int rc=0, lnam=0, lval=0;
   size_t i=0, n=0;

   if ( lst == (MetaParamList *) NULL || iobuf == (IO_BUFFER *) NULL )
      return -1;

   item_header.type = IO_TYPE_METAPARAM;           /* Metadata block is type 75 */

   if ( (rc = get_item_begin(iobuf,&item_header)) < 0 )
      return rc;

   if ( item_header.version != 0 )
   {
      fflush(stdout);
      fprintf(stderr,"Unsupported MetaParam version: %d.\n",
         item_header.version);
      get_item_end(iobuf,&item_header);
      return -1;
   }
   
   n = get_count(iobuf);
   
   clear_metaparam(lst);

   lst->ident = item_header.ident;
   lst->first = (MetaParamItem *) calloc(1,sizeof(MetaParamItem));
   item = lst->first;
   
   for (i=0; i<n && item != NULL; i++)
   {
      lnam = get_var_string(item_name,(int)sizeof(item_name),iobuf);
      lval = get_var_string(item_value,(int)sizeof(item_value),iobuf);
      if ( (size_t) lnam+1 >= sizeof(item_name) || (size_t) lval+1 >= sizeof(item_value) )
      {
         fflush(stdout);
         fprintf(stderr,"MetaParam entry at position %ju for ID=%ld might be truncated.\n", i+1, lst->ident);
      }
      item->name = strdup(item_name);
      item->value = strdup(item_value);
      if ( i+1 >= n )
         break;
      item->next = (MetaParamItem *) calloc(1,sizeof(MetaParamItem));
      item = item->next;
   }

   return get_item_end(iobuf,&item_header);
}

/* --------------------------- print_metaparam ------------------- */
/** @short Display the contents of a data block of meta parameters.
 *
 *  @param iobuf I/O buffer descriptor
 *
 *  @return 0 (OK), -1 or other non-zero (problem)
 */

int print_metaparam (IO_BUFFER *iobuf)
{
   IO_ITEM_HEADER item_header;
   char item_name[1024], item_value[10240];
   int rc=0;
   size_t i=0, n=0;

   if ( iobuf == (IO_BUFFER *) NULL )
      return -1;

   item_header.type = IO_TYPE_METAPARAM;           /* Metadata block is type 75 */

   if ( (rc = get_item_begin(iobuf,&item_header)) < 0 )
      return rc;

   if ( item_header.version != 0 )
   {
      fflush(stdout);
      fprintf(stderr,"Unsupported MetaParam version: %d.\n",
         item_header.version);
      get_item_end(iobuf,&item_header);
      return -1;
   }
   
   n = get_count(iobuf);
   
   printf("\nMeta parameters for ID=%ld%s with %ju name/value pairs.\n", 
      item_header.ident, (item_header.ident==-1)?" (global)":"", n);
   
   for (i=0; i<n; i++)
   {
      (void) get_var_string(item_name,sizeof(item_name)-1,iobuf);
      (void) get_var_string(item_value,sizeof(item_value)-1,iobuf);
      printf("   %s = %s\n", item_name, item_value);
      if ( i+1 >= n )
         break;
   }

   return get_item_end(iobuf,&item_header);
}

/* ---------------------- show_metaparam --------------------- */
/** @short Display the contents in a linked list of meta parameters,
 *         independent of its EventIO representation.
 *
 *  @param lst The starting point of a linked list. Must exist but can be without items.
 *
 *  @return 0 (OK), -1 or other non-zero (problem)
 */

int show_metaparam (const MetaParamList *lst)
{
   MetaParamItem *item = NULL, *next = NULL;

   if ( lst == (const MetaParamList *) NULL )
      return -1;

   printf("\nMeta parameters for ID=%ld.\n", lst->ident);

   /* Now we can write the items one by one */
   for ( item = lst->first; item != NULL; item = next)
   {
      next = item->next;
      if ( item->name == NULL )
         continue;
      printf("   %s = %s\n", item->name, item->value);
   }

   return 0;
}

/* ------------------------ search_metaparam -------------------- */
/** @short Search for a specific metaparameter by name (case-insensitive,
 *     first match returned) in a linked list of meta parameters.
 *     NULL is returned if there is no match.
 *
 *  Since metaparameters are stored only once per run (and telescope),
 *  and there should not be many of them, a linear search should be 
 *  good enough. 
 *
 *  @param lst The starting point of a linked list. Must exist but can be without items.
 *  @param name The name of the metaparameter searched for.
 *
 *  @return Value of metaparameter (first match found), or NULL pointer.
 */

const char *search_metaparam (const MetaParamList *lst, const char *name)
{
   MetaParamItem *item = NULL, *next = NULL;

   if ( lst == (const MetaParamList *) NULL || name == NULL )
      return NULL;

   for ( item = lst->first; item != NULL; item = next)
   {
      next = item->next;
      if ( item->name == NULL )
         continue;
      if ( strcasecmp(item->name,name) == 0 )
         return item->value;
   }

   return NULL;
}
