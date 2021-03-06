/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "smalloc.h"

#include "sparsematrix.h"




gmx_sparsematrix_t *
gmx_sparsematrix_init(int                    nrow)
{
    int i;
    gmx_sparsematrix_t *A;
    
    snew(A,1);
    
    A->nrow=nrow;
    snew(A->ndata,nrow);
    snew(A->nalloc,nrow);
    snew(A->data,nrow);
        
    for(i=0;i<nrow;i++) 
    {
        A->ndata[i]    = 0;
        A->nalloc[i]   = 0;
        A->data[i]     = NULL;
    }
    return A;    
}



void
gmx_sparsematrix_destroy(gmx_sparsematrix_t *   A)
{
    int i;
    
    /* Release each row */
    for(i=0;i<A->nrow;i++) 
    {
        if(A->data[i]!=NULL)
            sfree(A->data[i]);
    }
    /* Release the rowdata arrays */
    sfree(A->ndata);
    sfree(A->nalloc);
    sfree(A->data);
    /* Release matrix structure itself */
    sfree(A);    
}



void
gmx_sparsematrix_print(FILE *                 stream,
                       gmx_sparsematrix_t *   A)
{
    int i,j,k;
    
    for(i=0;i<A->nrow;i++)
    {
        if(A->ndata[i]==0) 
        {
            for(j=0;j<A->nrow;j++)
                fprintf(stream," %6.3f",0.0);
        }
        else
        {
            k=0;
            j=0;
            for(j=0;j<A->ndata[i];j++) 
            {
                while(k++<A->data[i][j].col) 
                    fprintf(stream," %6.3f",0.0);
                fprintf(stream," %6.3f",A->data[i][j].value);
            }
            while(k++<A->nrow)
                fprintf(stream," %6.3f",0.0);
        }
        fprintf(stream,"\n");
    }
    
}


real
gmx_sparsematrix_value(gmx_sparsematrix_t *    A,
                       int                     row, 
                       int                     col)
{
    bool found  = FALSE;
    int  i;
    real value;
    
    assert(row<A->nrow);

    value = 0;
    
    /* Find previous value */
    for(i=0;i<A->ndata[row] && (found==FALSE);i++)
    {
        if(A->data[row][i].col==col) 
        {
            found = TRUE;
            value = A->data[row][i].value;
        }
    }
    
    /* value=0 if we didn't find any match */
    return value;
}



void
gmx_sparsematrix_increment_value(gmx_sparsematrix_t *    A,
                                 int                     row, 
                                 int                     col,
                                 real                    difference)
{
    bool found  = FALSE;
    int i;
    
    assert(row<A->nrow);
    
    /* Try to find a previous entry with this row/col */
    for(i=0;i<A->ndata[row] && !found;i++)
    {
        if(A->data[row][i].col==col) 
        {
            found = TRUE;
            A->data[row][i].value += difference;
        }
    }
    
    /* Add a new entry if nothing was found */
    if(!found) 
    {
        /* add the value at the end of the row */
        if(A->ndata[row] == A->nalloc[row]) 
        {
            A->nalloc[row] += 100;
            if(A->data[row]==NULL)
            {
                snew(A->data[row],A->nalloc[row]);
            }
            else
            {
                srenew(A->data[row],A->nalloc[row]);
            }
        }
        A->data[row][A->ndata[row]].col   = col;
        /* Previous value was 0.0 */
        A->data[row][A->ndata[row]].value = difference;
        A->ndata[row]++;
    }     
}


/* Routine to compare column values of two entries, used for quicksort of each row.
 *
 * The data entries to compare are of the type gmx_sparsematrix_entry_t, but quicksort
 * uses void pointers as arguments, so we cast them back internally.
 */
static int
compare_columns(const void *v1, const void *v2)
{
    int c1 = ((gmx_sparsematrix_entry_t *)v1)->col;
    int c2 = ((gmx_sparsematrix_entry_t *)v2)->col;
    
    if(c1<c2)
        return -1;
    else if(c1>c2)
        return 1;
    else 
        return 0;
}


void
gmx_sparsematrix_compress(gmx_sparsematrix_t *    A)
{
    int i,j;
    
    for (i=0;i<A->nrow;i++) 
    {        
        /* Remove last value on this row while it is zero */
        while(A->ndata[i]>0 && A->data[i][A->ndata[i]-1].value==0)
            A->ndata[i]--;
        
        /* Go through values on this row and look for more zero elements */
        for(j=0;j<A->ndata[i];j++)
        {
            /* If this element was zero, exchange it with the last non-zero
             * element on the row (yes, this will invalidate the sort order)
             */
            if(A->data[i][j].value==0)
            {
                A->data[i][j].value = A->data[i][A->ndata[i]-1].value;
                A->data[i][j].col   = A->data[i][A->ndata[i]-1].col;
                A->ndata[i]--;
            }
        }
        /* Only non-zero elements remaining on this row. Sort them after column index */
        qsort((void *)(A->data[i]),
              A->ndata[i],
              sizeof(gmx_sparsematrix_entry_t),
              compare_columns);
    }        
}


void
gmx_sparsematrix_vector_multiply(gmx_sparsematrix_t *    A,
                                 real *                  x,
                                 real *                  y)
{
    real                        s,v,xi;
    int                         i,j,k;
    gmx_sparsematrix_entry_t *  data; /* pointer to simplify data access */
    
    for (i = 0; i < A->nrow; i ++)
        y[i] = 0;
    
    if(A->compressed_symmetric)
    {
        for (i = 0; i < A->nrow; i ++)
        {
            xi = x[i];
            s = 0.0;
            data = A->data[i];
            
            for (k=0;k<A->ndata[i];k++)
            {
                j = data[k].col;
                v = data[k].value;
                s += v * x[j];
                if(i!=j)
                    y[j] += v * xi; 
            }
            y[i] += s; 
        }    
    }
    else
    {
        /* not compressed symmetric storage */
        for (i = 0; i < A->nrow; i ++)
        {
            xi = x[i];
            s = 0.0;
            data = A->data[i];
            
            for (k=0;k<A->ndata[i];k++) 
            {
                j = data[k].col;
                v = data[k].value;
                s += v * x[j];
            }
            y[i] += s; 
        }    
    }
}



