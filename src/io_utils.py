## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## 
## 	Utilities for reading and writing files for kPadre
## 
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def time_report():
	'''Returns date and time'''

	from datetime import datetime
	now = datetime.now()
	current_time = now.strftime("[%y/%m/%d %H:%M:%S]: ")
	return current_time

def write_scores_table(assignments, fname):
	'''Writes a text file with the results of the analysis of
	   the AssignGenome() function '''
 
	with open(f'KP_{fname}_readAssignment.txt', 'w') as fh_out:
 
		fh_out.write('ReadID\tmaxG1\tmaxG2\tgenome\tchr_G1\tpos_G1\tchrG2\tposG2\tgroup\n')

		for k,v in assignments.items():
			read = k
			genome = v[-1][0]
			maxG1  = v[-1][1]
			maxG2  = v[-1][2]
			chrG1  = v[-1][3]
			posG1  = v[-1][4]
			chrG2  = v[-1][5]
			posG2  = v[-1][6]
			group  = v[-1][7]
	 
			fh_out.write(f'{read}\t{maxG1}\t{maxG2}\t{genome}\t{chrG1}\t{posG1}\t{chrG2}\t{posG2}\t{group}\n')
