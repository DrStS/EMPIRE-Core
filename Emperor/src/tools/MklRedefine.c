/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Munich
 *
 *  All rights reserved.
 *
 *  This file is part of EMPIRE.
 *
 *  EMPIRE is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EMPIRE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with EMPIRE.  If not, see http://www.gnu.org/licenses/.
 */
/*#include <stdio.h>
#include <string.h>
#define BUFLEN 16
int mkl_progress_( int* ithr, int* step, char* stage, int lstage )
{
  char buf[BUFLEN];
  if( lstage >= BUFLEN ) lstage = BUFLEN-1;
  strncpy( buf, stage, lstage );
  buf[lstage] = '\0';
  printf( "In thread %i, at stage %s, steps passed %i\n", *ithr, buf, *step );
  return 0;
}*/
