target= "16S ribosomal RNA"
results = []
result_file = open("all_16s.fasta", "w+")

def one(name,text):
	for item in text.split(">"):
			if target in item:
					print(name)
					text = ">"+item.replace(target,name + " " + target)
					result_file.write(text)
for file in open("target.txt","r").read().split("\n"):
	path = "/working_dir/work/prokka_result/"+file+"/PROKKA_11162023.ffn"
	try:
		file2 = open(path,"r").read()
		one(file, file2)
	except:
		pass