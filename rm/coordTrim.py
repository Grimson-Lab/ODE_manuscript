"""-----------------------------------------------------------------------------
  * Script Name: coordTrim.py
  * Description: This script will iterate through the output of any .bed file that is
	an output of your previous "controlMerge.py" scripts or that has a .bed format 
	alongside a divergence value in column 8. Iterate through and shrink TE coords 
	based on the divergence score (lower wins) across all intersections. Not as 
	fast as it could be, but considering it handles ~ 4k intervals a second, you've
	written worse code
  * Created By:  Jason David Chobirko
  * Date:        Thursday December 25th, 2025
-----------------------------------------------------------------------------"""
# First import all of the necessary modules you need
import sys, os, re
import pandas as pd

# The expected call to run the script is:
# python coordTrim.py input.bed min.size output.bed

# Use pandas so that you can iterate through and change each row/column as needed. Slow :'(
rf = pd.read_csv(sys.argv[1], sep='\t', header=None, names=['chr', 'start', 'end', 'name', 'comb', 'strand', 'con', 'div'])

# Now to iterate through each row in the above object
for x in range(0,len(rf)-1):
	# Ignore entries that have a current length of less than provided value, defaulted to 50bps
	if (rf.loc[x,'end'] - rf.loc[x,'start'] < int(sys.argv[2])):
	# if (rf.loc[x,'end'] - rf.loc[x,'start'] < 50):
		continue
	# Check to see if the current entry overlaps the next entry. Keep going until there are no more intersects
	count = 1
	while (rf.loc[x,'chr'] == rf.loc[x+count,'chr'] and rf.loc[x,'end'] > rf.loc[x+count,'start']):
		# The next element is intersecting. If the current element is more diverged, shrink. Otherwise, shrink the next element
		# up to its end coordinates and increment count to check the next element in the list. Make sure count doesn't go over the 
		# full length of rf
		if (rf.loc[x,'div'] <= rf.loc[x+count,'div']):
			rf.loc[x+count,'start'] = min(rf.loc[x+count,'end'],rf.loc[x,'end'])
			# rf.loc[x+count,'edit'] += 1
			count += 1
		else:
			# In the *RARE* case that a full-length element fully encompasses a smaller less-diverged element (such as ORR1B1-int 
			# and a SINE element, for example, remove the smaller element instead
			if (rf.loc[x,'end'] - rf.loc[x,'start'] >= 2.5 * (rf.loc[x+count,'end'] - rf.loc[x+count,'start']) and rf.loc[x,'end'] - rf.loc[x,'start'] >= 400):
				rf.loc[x+count,'start'] = min(rf.loc[x+count,'end'],rf.loc[x,'end'])
				# rf.loc[x+count,'edit'] += 1
				count += 1
			else:
				rf.loc[x,'end'] = min(rf.loc[x+count,'start'],rf.loc[x,'end'])
				# rf.loc[x,'edit'] += 1
		# Check to make sure that count + current row don't go beyond the bounds of the data.frame
		if (x+count == len(rf)):
			break

# Finally, remove all the entries that are less than or equal to the input size and write to file!
rf = rf[rf['end'] - rf['start'] > int(sys.argv[2])]
rf.sort_values(['chr', 'start'], ascending = [True, True])
rf.to_csv(sys.argv[3], sep ='\t', mode='a', index=False, header=False)

