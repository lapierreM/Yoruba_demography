#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
	Copyright (C) 2016  Lapierre/Marguerite

	This program is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public License
	as published by the Free Software Foundation; either version 2.1
	of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

	for more information, please contact: marguerite.lapierre@mnhn.fr
"""

#------------------------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------------------------

# Open file in which the SFS will be written in the format compatible with dadi
f_dadi=open("YRI.fs","w")

# Read the original SFS file
with open("../../data/FoldSFS_YRI.txt","r") as f_sfs:
	f_dadi.write("217 folded \"YRI\"\n") # dadi header
	f_dadi.write("0 ") # value for frequency=0 (does not matter)
	
	for ligne in f_sfs:
		liste=ligne.split('\t')
		f_dadi.write(liste[1][:-1]+" ") # folded SFS values

	for i in xrange(107):
		f_dadi.write("0 ")  # values for frequency > n/2 (folded SFS so 0)

f_dadi.write("0\n1 ") # value for frequency=1 (does not matter), newline, 1 to mask the value for frequency=0

for i in xrange(108):
	f_dadi.write("0 ") # 0 to unmask the values of the folded SFS
for i in xrange(107):
	f_dadi.write("1 ") # 1 to mask the values for "unfolded" SFS
f_dadi.write("1\n") # 1 to mask the value for frequency=1

f_dadi.close()