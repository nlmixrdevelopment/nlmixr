ss = [\
["CL", "V", "KA", "TLAG"],
["CL", "V", "CLD", "VT", "KA", "TLAG"],
["CL", "V", "CLD", "VT", "CLD2", "VT2", "KA", "TLAG"],
["KE", "V", "KA", "TLAG"],
["KE", "V", "K12", "K21", "KA", "TLAG"],
["KE", "V", "K12", "K21", "K13", "K31", "KA", "TLAG"]
]

for i in ss:
	s = i[:-2]
	s1 = ['l'+x for x in s]
	s2 = [x+' = exp(l%s)'%x for x in s]
	print 'function(%s){\nc(\n%s\n)\n}\n' %(', '.join(s1), ',\n'.join(s2))
	s = i[:-1]
	s1 = ['l'+x for x in s]
	s2 = [x+' = exp(l%s)'%x for x in s]
	print 'function(%s){\nc(\n%s\n)\n}\n' %(', '.join(s1), ',\n'.join(s2))
	s = i
	s1 = ['l'+x for x in s]
	s2 = [x+' = exp(l%s)'%x for x in s]
	print 'function(%s){\nc(\n%s\n)\n}\n' %(', '.join(s1), ',\n'.join(s2))
