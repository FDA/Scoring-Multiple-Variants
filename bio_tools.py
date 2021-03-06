import tools
import os, re, sys, requests, time, copy, bio_tools, warnings, urllib.parse, traceback, warnings, copy, xlrd
from bs4 import BeautifulSoup
from Bio import Entrez
import parse_range
import pandas as pd
import numpy as np
import score_variants, compute_features
Entrez.email=""

complement = {"A":"T", "T":"A", "G":"C", "C":"G"}

chromosome_map = {"GRCh37":{"1":"NC_000001.10", "2":"NC_000002.11", "3":"NC_000003.11", "4":"NC_000004.11", "5":"NC_000005.9", "6":"NC_000006.11", "7":"NC_000007.13", "8":"NC_000008.10", "9":"NC_000009.11", "10":"NC_000010.10", "11":"NC_000011.9", "12":"NC_000012.11", "13":"NC_000013.10", "14":"NC_000014.8", "15":"NC_000015.9", "16":"NC_000016.9", "17":"NC_000017.10", "18":"NC_000018.9", "19":"NC_000019.9", "20":"NC_000020.10", "21":"NC_000021.8", "22":"NC_000022.10", "X":"NC_000023.10", "Y":"NC_000024.9"},"GRCh38":{"1":"NC_000001.11", "2":"NC_000002.12", "3":"NC_000003.12", "4":"NC_000004.12", "5":"NC_000005.10", "6":"NC_000006.12", "7":"NC_000007.14", "8":"NC_000008.11", "9":"NC_000009.12", "10":"NC_000010.11", "11":"NC_000011.10", "12":"NC_000012.12", "13":"NC_000013.11", "14":"NC_000014.9", "15":"NC_000015.10", "16":"NC_000016.10", "17":"NC_000017.11", "18":"NC_000018.10", "19":"NC_000019.10", "20":"NC_000020.11", "21":"NC_000021.9", "22":"NC_000022.11", "X":"NC_000023.11", "Y":"NC_000024.10"}}

#if sequences names are given as NM_######:###-###, will merge separated subsequences for alignment
def merge_seqs(seqs):
	seqs = {k.split()[0]:v for k,v in seqs.items()}
	seqs = {k.split(":")[0]:{k2.split(":")[1]:v2 for k2,v2 in seqs.items() if k2.split(":")[0] == k.split(":")[0]} for k,v in seqs.items()}
	seqs = {k:sorted(v.items(), key=lambda kv: int(kv[0].split("-")[0])) for k,v in seqs.items()}
	seqs = {k:"".join([kv[1] for kv in v]) for k,v in seqs.items()}
	return(seqs)

#returns a refseq/genbank accession without the version
def acc_wo_version(acc):
	return(".".join(acc.split(".")[:-1]))

#determines if an accession is contained in a list of accessions, ignoring version number
def compare_acc_list(acc, l):
	return(acc_wo_version(acc) in [acc_wo_version(x) for x in l])

def find_acc_match(acc, l):
	l = [li for li in l if acc_wo_version(li) == acc_wo_version(acc)]
	l = sorted(l, reverse=True, key=lambda k: int(k.split(".")[-1]))
	return(l[0])

#find the chromosome and location of a gene (using standard gene symbol)
def get_genome_loc(genename, assembly="GRCh37", organism="homo sapiens"):
	r = tools.retry_request(requests.get, ["https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term="+genename+"[gene]+AND+"+re.sub("\s+", "+", organism, flags=re.MULTILINE)+"[organism]"])
	ids = re.findall("<Id>(\d+)</Id>", r.text)
	for accid in ids:
		try:
			handle = tools.retry_func(Entrez.efetch, [], {"db":"gene", "id":accid, "retmode":"xml"})
			record = handle.read()
			if genename in str(record):
				chromosome = max(re.findall("Chromosome (\w+)", str(record)), key=len)
				gid = chromosome_map[assembly][str(chromosome)]
				spl = [x for x in record.split("<Gene-commentary>") if gid.split(".")[0] in x]
				spl = [x for x in spl if re.search("\<Gene\-commentary\_version\>" + gid.split(".")[-1] + "\<\/Gene\-commentary\_version\>", x)]
				if len(spl) > 0:
					s = spl[0]
					r1 = int(max(re.findall("\<Seq\-interval\_from\>(\d+)\<\/Seq\-interval\_from\>", s), key=len))
					r2 = int(max(re.findall("\<Seq\-interval\_to\>(\d+)\<\/Seq\-interval\_to\>", s), key=len))
					return(chromosome, [r1, r2])
		except Exception as e:
			print(e)
	return(None)

#get all accession ids for a gene name (transcript, genomic, etc)
def get_accids(genename, organism="homo sapiens"):
	r = tools.retry_request(requests.get, ["https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term="+genename+"[gene]+AND+"+re.sub("\s+", "+", organism, flags=re.MULTILINE)+"[organism]"])
	ids = re.findall("<Id>(\d+)</Id>", r.text)
	for accid in ids:
		try:
			handle = tools.retry_func(Entrez.efetch, [], {"db":"gene", "id":accid, "retmode":"xml"})
			record = str(handle.read().decode("utf-8"))
			if genename in str(record):
				accids = re.findall("\<Gene\-commentary\_accession\>(\S+)\<\/Gene\-commentary\_accession\>\s+\<Gene\-commentary\_version\>(\d+)\<\/Gene\-commentary\_version\>", record)
				accids = [li[0] + "." + li[1] for li in accids]
				genomic = [li for li in accids if li.startswith("NG_")]
				mRNA = [li for li in accids if li.startswith("NM_")]
				protein = [li for li in accids if li.startswith("NP_")]
				if len(mRNA) != len(protein):
					raise Warning("mRNA and protein accession ids may not line up")
				mRNA = [(mRNA[i], protein[i]) for i in range(len(mRNA))]
				if len(genomic) > 0 and len(mRNA) > 0:
					genomic = sorted(list(set(genomic)), key=lambda kv: int(max(re.findall("\d+", kv), key=len)))
					mRNA = sorted(list(set(mRNA)), key=lambda kv: int(max(re.findall("\d+", kv[0]), key=len)))
					return({"genomic":genomic, "mRNA":mRNA})
		except Exception as e:
			print(e)
	try:
		return({"genomic":genomic, "mRNA":mRNA})
	except Exception as e:
		print(e)

#get the sequences corresponding to a gene, given the genomic and transcript accession ids
def get_seqs_and_exons(genename, genomic, mRNA, organism="homo sapiens", index=1):
	transcripts = {}
	mRNAs = {}
	for genid in genomic:
		try:
			handle = tools.retry_func(Entrez.efetch, [], {"db":"nucleotide", "id":genid, "retmode":"xml"})
			genomic_record = str(handle.read().decode("utf-8"))
			genomic_seq = max(re.findall("\<GBSeq\_sequence\>([acgtACGT\s]+)\<\/GBSeq\_sequence\>", genomic_record), key=len)
			genomic_seq = re.sub("\s+", "", genomic_seq).upper()
			for piece in genomic_record.split("<GBFeature>"):
				if "<GBFeature_key>mRNA</GBFeature_key>" in piece:
					transcript_id = max(re.findall("\<GBQualifier\_value\>(NM\_.+)\<\/GBQualifier\_value\>", piece), key=len)
					transcripts[transcript_id] = max(re.findall("\<GBFeature\_location\>(.+)\<\/GBFeature\_location\>", piece), key=len)
				elif "<GBFeature_key>CDS</GBFeature_key>" in piece:
					mRNA_id = max(re.findall("\<GBQualifier\_value\>(NP\_.+)\<\/GBQualifier\_value\>", piece), key=len)
					mRNAs[mRNA_id] = max(re.findall("\<GBFeature\_location\>(.+)\<\/GBFeature\_location\>", piece), key=len)
		except Exception as e:
			print(e)
			traceback.print_exc()
	if len(genomic) < 1: #no genomic ids found
		for i in range(len(mRNA)):
			transcripts[mRNA[i][0]] = None
			mRNAs[mRNA[i][1]] = None
		genomic_seq = None
	for (mRNA_acc, protein_acc) in mRNA:
		if compare_acc_list(mRNA_acc,transcripts) and compare_acc_list(protein_acc,mRNAs):
			try:
				time.sleep(1)
				handle = tools.retry_func(Entrez.efetch, [], {"db":"nucleotide", "id":mRNA_acc, "retmode":"xml"})
				record = str(handle.read().decode("utf-8"))
				seq = max(re.findall("\<GBSeq\_sequence\>([acgtACGT\s]+)\<\/GBSeq\_sequence\>", record), key=len)
				seq = re.sub("\s+", "", seq).upper()
				for piece in record.split("<GBFeature>"):	
					if re.search("\<GBQualifier\_value\>NP\_.+\<\/GBQualifier\_value\>", piece):
						cds = max(re.findall("\<GBFeature\_location\>(.+)\<\/GBFeature\_location\>", piece), key=len)
						ends = [int(x) for x in re.findall("\d+", cds)]
						if len(ends) == 2:
							seqs = {"genomic":genomic_seq, "transcript":seq, "ORF":seq[ends[0]-index:ends[1]]}
							cdses = {"ORF":mRNAs[find_acc_match(protein_acc, mRNAs)], "transcript":transcripts[find_acc_match(mRNA_acc, transcripts)]}
							return(seqs, cdses, mRNA_acc, protein_acc)
					
			except Exception as e:
				print(e)
				traceback.print_exc()

#find the genomic coordinates of a transcript coordinate mutation
def cds_to_genomic(accid, ntvar, assembly="GRCh37"):
	version = int(accid.split(".")[-1])
	accid_wo_version = ".".join(accid.split(".")[:-1])
	while version >= 0: #mutalyzer may not have up-to-date versions, so keep looking until you find one
		try:
			r = requests.get("https://mutalyzer.nl/position-converter?assembly_name_or_alias="+assembly+"&description="+accid_wo_version+"."+str(version)+"%3Ac."+str(ntvar[0])+ntvar[1][0]+"%3E"+ntvar[1][1])
			genvar = re.findall("\<h4\>Chromosomal variant\<\/h4\>\s*\<p\>\<code\>(.+)\<\/code\>\<\/p\>", r.text)
		except Exception as e:
			version -= 1
		if genvar == None or len(genvar) == 0:
			version -= 1
		else:
			break
	if genvar != None and len(genvar) > 0:
		genvar = [x.replace("&gt;", ">") for x in genvar]
		return(genvar[0])
	else:
		return(None)
def rev_chromosome_map(nc, assembly="GRCh37"):
	k = [k for k,v in chromosome_map[assembly].items() if v==nc]
	if len(k) > 0:
		k = k[0]
		if tools.is_numeric(k):
			return(k)
		elif k == "X":
			return('23')
		else:
			return(None)
	else:
		return(None)

def find_freqs(ids, gene=""):
	data = {}
	init_time = time.time()
	for i, idi in enumerate(ids):
		server = "https://rest.ensembl.org"
		try:
			ext = "/vep/human/id/" + idi
			r = tools.retry_request(requests.get, positional_arguments=[server+ext], keyword_arguments={"headers":{"Content-Type":"application/json"}})
			decoded = r.json()
			for j in range(len(decoded[0]["transcript_consequences"])):
				sub = decoded[0]["transcript_consequences"][j]
				print(sub)
				if "gene symbol" in sub and sub["gene_symbol"] == gene:
					if "cds_start" in sub and "cds_end" in sub and sub["cds_start"] == sub["cds_end"]:
						pos = int(sub["cds_start"])
						break
			for j in range(len(decoded[0]["colocated_variants"])):
				try:
					data[pos] = decoded[0]["colocated_variants"][j]["frequencies"]
					break
				except KeyError as e:
					pass
				except Exception as e:
					print(e)
					traceback.print_exc()
		except Exception as e:
			print(e)
			traceback.print_exc()
		tools.update_time(i, len(ids), init_time)
	return(data)

def convert_rsid_to_pos(ids, refid):
	data = {}
	refid = refid.replace("_", "\_").replace(".", "\.")
	for rsid in ids:
		rsid = max(re.findall("\d+", rsid), key=len)
		s = tools.retry_func(Entrez.efetch, [], {'db':'snp', 'id':rsid, 'retmode':'xml'}).read()
		try:
			loc = re.findall(refid + "\.?\d*\:g\.\d+[ACGT]\&gt\;[ACGT]",s)
			data[rsid] = loc[0]
		except IndexError as e:
			pass
		except Exception as e:
			print(e)
			traceback.print_exc()
		finally:
			time.sleep(2)
	return(data)
	
def cross_ref(df1, df2, refid="NM_024006.6"):
	data = {}
	for i, row in df1.iterrows():
		ids = tools.get_mutalyzer(refid, i[:-1] + ">" + i[-1], "c.")
		for i in ids:
			refid = i.split(":")[0]
			pos = int(max(re.findall("\d+", i.split(":")[-1]), key=len))
			if refid == "NC_000016.9":
				print(df2.loc[df2["POS"] == pos].T.to_dict())
				data.update(df2.loc[df2["POS"] == pos].T.to_dict())
	return(data)
	
def transcript_to_assembly(var, assembly="GRCh37"):
	server = "https://rest.ensembl.org"
	ext = "/variant_recoder/human/" + var
	r = tools.retry_request(requests.get, [server+ext], {'headers':{ "Content-Type" : "application/json"}})
	accids = re.findall("\"([\S]+?)\"", r.text)
	accids = [accid for accid in accids if accid.startswith("NC")]
	accids = [accid for accid in accids if accid.count(":") < 3]
	positions = [int(max(re.findall("\d+", "".join(accid.split(":")[1:])), key=len)) for accid in accids]
	nc = [accid.split(":")[0] for accid in accids]
	if len(positions) != len(nc):
		raise Exception("Problem parsing variant positions " + "\t" + str(accids))
	for i in range(len(positions)):
		for k,v in chromosome_map.items():
			if nc[i] in v.values():
				ch = [k2 for k2 in v.keys() if v[k2] == nc[i]][0]
				if k == assembly:
					return((accids[i], nc[i], positions[i], ch))
				else:
					ext = "/map/human/" + k + "/" + str(ch) + ":" + str(positions[i]) + ".." + str(positions[i]) + "/" + assembly + "?"
					r = tools.retry_request(requests.get, [server+ext], {'headers':{ "Content-Type" : "application/json"}})
					r = r.json()
					for i in range(len(r["mappings"])):
						try:
							pos = r["mappings"][i]["mapped"]['start']
							nc = chromosome_map[assembly][str(ch)]
							nts = re.findall("[A-Za-z]", var)
							accid = str(nc) + ":g." + str(pos) + nts[-2] + ">" + nts[-1]
							return((accid, nc, pos, ch))
						except Exception as e:
							print(e)
							traceback.print_exc()
	return((None, None, None, None))

def convert_assembly(var, assembly="GRCh38"):
	server = "https://rest.ensembl.org"
	nc = var.split(":")[0]
	position = int(max(re.findall("\d+", "".join(var.split(":")[1:])), key=len))
	for k,v in chromosome_map.items():
		if nc in v.values():
			ch = [k2 for k2 in v.keys() if v[k2] == nc][0]
			if k == assembly:
				return((var, nc, position, ch))
			else:
				ext = "/map/human/" + k + "/" + str(ch) + ":" + str(position) + ".." + str(position) + "/" + assembly + "?"
				r = tools.retry_request(requests.get, [server+ext], {'headers':{ "Content-Type" : "application/json"}})
				r = r.json()
				for i in range(len(r["mappings"])):
					try:
						pos = r["mappings"][i]["mapped"]['start']
						nc = chromosome_map[assembly][str(ch)]
						nts = re.findall("[A-Za-z]", var)
						accid = str(nc) + ":g." + str(pos) + nts[-2] + ">" + nts[-1]
						return((accid, nc, pos, ch))
					except Exception as e:
						print(e)
						traceback.print_exc()
	return((None, None, None, None))



def transcript_to_gene_info(var, email=Entrez.email):
	Entrez.email = email
	handle = tools.retry_func(Entrez.efetch, [], {'db':"nucleotide", 'id':var.split(":")[0], 'rettype':"gb", 'retmode':"text"})
	record = handle.read()
	try:
		gene = max(re.findall("\/gene\=\"(DGKB)\"\\n", record), key=len)
	except Exception as e:
		gene = None
	try:
		uniprot = max(re.findall("UniProtKB\:[A-Z|a-z|\d]+", record), key=len)
	except Exception as e:
		uniprot = None
	return(gene, uniprot)

def assembly_to_transcript(var, nm=None):
	server = "https://rest.ensembl.org"
	ext = "/variant_recoder/human/" + var
	r = tools.retry_request(requests.get, [server+ext], {'headers':{ "Content-Type" : "application/json"}})
	accids = re.findall("\"([\S]+?)\"", r.text)
	accids = [accid for accid in accids if accid.startswith("NM")]
	if nm == None:
		nm = min(set([accid.split(":")[0] for accid in accids]), key=lambda k: int(max(re.findall("\d+", k.split(".")[0]), key=len)))
	accids = [accid for accid in accids if nm in accid]
	return(accids)	

def get_vars(df, ch, r, studyname = "", index_col=None, chromosome="#CHR", position="POS"):
	if isinstance(df, str):
		if studyname == "" or studyname == None:
			studyname = df.split("_")[2]
		if index_col != None:
			df = pd.read_csv(df, sep="\t", index_col=index_col)
		else:
			df = pd.read_csv(df, sep="\t")
	data = {}
	if tools.is_numeric(ch):
		df = df.loc[df[chromosome] == int(ch)]
	else:
		if ch.upper() == "X":
			df = df.loc[df[chromosome] == 23]
		if ch.upper() == "Y":
			return({})
	df = df.loc[df[position] >= r[0]-6000]
	df = df.loc[df[position] <= r[1]+6000]
	df = copy.deepcopy(df.T.to_dict())
	try:
		for k,v in df.items():
			df[k]["Study"] = studyname
	except Exception as e:
		print(e)
	return(df)

def batch_mutalyzer(df, mutalyzer, assembly="GRCh37"):
	trans = {}
	with open(mutalyzer, "r") as inf:
		lines = inf.readlines()
		for line in lines:
			line = [line.split("\t")[0], "\t".join(line.split("\t")[1:])]
			line[1] = re.findall("[A-Za-z]+\_\d+\.?\d*\:[a-z]\.[\-|\+]?\d+[\-|\+]?\d*[A-Za-z]\&gt\;[A-Za-z]", line[1])+re.findall("[A-Za-z]+\_\d+\.?\d*\:[a-z]\.\-?\d+\-?\d*[A-Za-z]\>[A-Za-z]", line[1])
			trans[line[0]] = line[1]
	df = pd.read_csv(df, sep="\t", index_col=4)
	d = {}
	for k,v in trans.items():
		try:
			nc = k.split(":")[0]
			ch = rev_chromosome_map(nc, assembly)
			var = k.split(".")[-1]
			pos = tools.get_num(var)
			alleles = re.findall("[A-Z]+", var)
			tmpdf = df.loc[[str(ch) + ":" + str(pos) + ":" + str(alleles[0]) + ":" + str(alleles[1])]]
			tmpdf = tmpdf.T.to_dict()
			for k in tmpdf.keys():
				tmpdf[k]["Trans"] = [x for x in v if x.startswith("NM")]
			d.update(tmpdf)
		except Exception as e:
			print(k + "\t" + str(e))
			traceback.print_exc()
	return(d)

def get_genes(d):
	refids = []
	for k in d.keys():
		if isinstance(d[k]["Trans"], (list, tuple,)):
			for x in d[k]["Trans"]:
				refids.append(x.split(":")[0])
	refids = {x:None for x in refids}
	init_time = time.time()
	for i, k in enumerate(refids):
		handle = tools.retry_func(Entrez.efetch, [], {'db':'nucleotide', 'id':k, 'retmode':'xml'})
		record = handle.read()
		names = re.findall("\/gene=\"(\S+)\"", record)
		if len(names) < 1:
			names = re.findall("\<GBQualifier\_name\>gene\<\/GBQualifier\_name\>\s+\<GBQualifier\_value\>(.+)\<\/GBQualifier\_value\>", record)
		names = list(set(names))
		refids[k] = names
		time.sleep(2)
		tools.update_time(i, len(refids), init_time)
	for k in d.keys():
		try:
			d[k]["GENE"] = [refids[refid] for refid in [x.split(":")[0] for x in d[k]["Trans"]] if refid in refids]
			d[k]["GENE"] = list(set(set().union(*d[k]["GENE"])))
			if len(d[k]["GENE"]) == 1:
				d[k]["GENE"] = d[k]["GENE"][0]
			elif len(d[k]["GENE"]) < 1:
				d[k]["GENE"] = None
		except Exception as e:
			print(e)
			traceback.print_exc()
	return(d)

def gwas_pipeline(dffile, genelist, prefix=None, index_col="SNP", p_col="all_inv_var_meta_p", chromosome="#CHR", position="POS", ref="REF", alt="ALT", plim=0.05, translation=None, assembly="GRCh37", check=False):
	nm_ids = []
	for gene in genelist:
		accids = bio_tools.get_accids(gene)
		nm_ids += [xi[0] for xi in accids["mRNA"]]
	nm_ids = list(set(nm_ids))
	if prefix == None:
		prefix = os.path.splitext(dffile)[0]
	df = pd.read_csv(dffile, sep="\t", index_col=index_col)
	df = df.loc[df[p_col] <= plim]
	newdf = {}
	for gene in genelist:
		try:
			loc = bio_tools.get_genome_loc(gene, assembly=assembly)
			newdf.update(get_vars(df, loc[0], loc[1], chromosome=chromosome, position=position))
		except TypeError as e:
			pass
	with open(prefix + "_nc_acc.txt", "w") as outf:
		for k,v in newdf.items():
			if v[chromosome] < 23:
				nc = chromosome_map[assembly][str(v[chromosome])]
			elif v[chromosome] == 23:
				nc = chromosome_map[assembly]['X']
			elif v[chromosome] == 23:
				nc = chromosome_map[assembly]['Y']
			newdf[k]["Genomic"] = nc + ":g." + str(v[position]) + v[ref] + ">" + v[alt]
			outf.write(newdf[k]["Genomic"] + "\n")
	if check:
		with open(prefix + "_nc_acc.txt", "w") as outf:
			for k,v in list(newdf.items()):
				var = v["Genomic"].split(".")[-1]
				nc_id = v["Genomic"].split(":")[0]
				errors, warnings = tools.check_variant(nc_id, var, "g.")
				if max(errors) > 0:
					nts = re.findall("[A-Z]+", var)
					pos = str(tools.get_num(var))
					newerrors, newwarnings = tools.check_variant(nc_id, pos + nts[1] + ">" + nts[0], "g.")
					if max(newerrors) < 1: #flip because study mixes up WT and mut alleles (inconsistent usage)
						newvar = nc_id + ":g." + pos + nts[1] + ">" + nts[0]
						newname = "chr" + [k for k in chromosome_map[assembly].keys() if chromosome_map[assembly][k] == nc_id][0] + ":" + pos + ":" + nts[1] + ":" + nts[0]
						newdf[newname] = copy.deepcopy(v)
						newdf[newname][ref] = nts[1]
						newdf[newname][alt] = nts[0]
						newdf[newname]["Genomic"] = newvar
						outf.write(newvar + "\t")
						del newdf[k]
					else:
						warnings.warn("Mismatch " + str(k) + "\t" + str(v["Genomic"]))
				else:
					outf.write(v["Genomic"] + "\n")
	trans = {}
	if translation != None and not os.path.exists(translation):
		input("Press Enter to continue...")
	try:
		with open(translation, "r") as outf:
			for i, line in enumerate(outf.readlines()):
				if i > 0:
					trans[line.split("\t")[0]] = line.strip().split("\t")[1:]
	except:
		for k,v in newdf.items():
			trans[v["Genomic"]] = get_mutalyzer(v["Genomic"].split(":")[0], v["Genomic"].split(".")[-1], "g.", assembly=assembly)
	trans = {k:[x for x in v if x.startswith("NM")] for k,v in trans.items()}
	#trans = {k:[x for x in v if bio_tools.compare_acc_list(x.split(":")[0], nm_ids)] for k,v in trans.items()}
	trans = {k:v for k,v in trans.items() if len(v) > 0}
	trans = {k:min(v, key=lambda kv: tools.get_num(kv.split(".")[0])) for k,v in trans.items()}
	print(trans)
	for (k,v) in list(newdf.items()):
		if v["Genomic"] not in trans:
			del newdf[k]
	seqs = get_subseq(trans)
	for k,v in newdf.items():
		try:
			newdf[k]["Transcript"] = trans[v["Genomic"]]
			newdf[k].update(seqs[trans[v["Genomic"]]])
			newdf[k].update(compare_splicing(seqs[trans[v["Genomic"]]]["WT"], seqs[trans[v["Genomic"]]]["MUT"]))
		except Exception as e:
			print(e)
			traceback.print_exc()
	for k,v in newdf.items():
		newdf[k].update(get_all_variant_info(v["Genomic"]))
	pd.DataFrame(newdf).T.to_csv(prefix + "_genes.tsv", sep="\t")
	return(newdf)

def get_subseq(trans):
	trans = {k:v for k,v in trans.items() if v != None and (isinstance(v, (list, tuple,)) or isinstance(v, str))}
	for k,v in trans.items():
		if isinstance(v, (list, tuple,)) and not isinstance(v, str):
			try:
				trans[k] = [vi for vi in v if vi.startswith("NM")][0]
			except IndexError as e:
				trans[k] = v[0]
	revtrans = {v:k for k,v in trans.items()}
	nm_accs = set([nmvar.split(".")[0] for nmvar in trans.values()])
	nm_accs = {k:[v for v in trans.values() if v.split(".")[0] == k] for k in nm_accs}
	muts = {}
	for k,v in nm_accs.items():
		gene = get_genes({v[0]:{"Trans":[v[0]]}})[v[0]]['GENE']
		if isinstance(gene, (list, tuple,)) and len(gene) > 1:
			gene = gene[0]
		elif gene == None or len(gene) < 1:
			warnings.warn("No gene information found")

		r = [[np.inf, 0, ""], [-np.inf, 0, ""]]
		for nm_var in v:
			position = max(re.findall("[\-|\+]?\d+[\-|\+]?\d*", nm_var.split(".")[-1]), key=len)
			position = re.findall("[\-|\+]?\d+", position)
			if len(position) == 1:
				position = [int(position[0]), 0]
			else:
				position = [int(position[0]), int(position[1])]
			if position[0] < r[0][0] or (position[0] == r[0][0] and position[1] < r[0][1]):
				r[0] = [position[0], position[1], nm_var]
			if position[0] > r[1][0] or (position[0] == r[1][0] and position[1] > r[1][1]):
				r[1] = [position[0], position[1], nm_var]
		try:
			seqs, cdses, mRNA, protein = bio_tools.pipeline(gene)
		except Exception as e:
			continue
		flip = False
		minpos = revtrans[r[0][2]]
		maxpos = revtrans[r[1][2]]
		ncmin = minpos.split(":")[0]
		ncmax = maxpos.split(":")[0]
		if ncmin != ncmax:
			warnings.warn("NC versions don't match: " + str(ncmin) + "\t" + str(ncmax))
		minpos2 = min([minpos, maxpos], key=lambda kv: tools.get_num(kv.split(".")[-1]))
		maxpos2 = max([minpos, maxpos], key=lambda kv: tools.get_num(kv.split(".")[-1]))
		if maxpos2 == minpos and minpos2 == maxpos:
			flip = True
		minpos = tools.get_num(minpos2.split(":")[-1])
		maxpos = tools.get_num(maxpos2.split(":")[-1])
		handle = tools.retry_func(Entrez.efetch, [], {'db':"nucleotide", 'id':ncmin, 'rettype':"fasta", 'retmode':"text", 'seq_start':minpos-250, 'seq_stop':maxpos+250})
		record = handle.read()
		seq = re.sub("\s+", "", "".join(record.split("\n")[1:]))
		for z, nm_var in enumerate(v):
			try:
				nc_var = revtrans[nm_var]
				nm_change = max(re.findall("\d+([A-Z]+\>[A-Z]+)", nm_var), key=len)
				nc_change = max(re.findall("\d+([A-Z]+\>[A-Z]+)", nc_var), key=len)
				nm_change = re.findall("[A-Z]+", nm_change)
				nc_change = re.findall("[A-Z]+", nc_change)
				minus = []
				mixup = False
				for i in range(len(nm_change)):
					if nm_change[i] == nc_change[i]:
						minus.append(False)
					elif "".join([complement[nm_change[i][j]] for j in range(len(nm_change[i]))]) == nc_change[i]:
						minus.append(True)
					else:
						minus.append(False)
						mixup=True
				minus = all(minus)
				if mixup:
					if nc_var not in get_mutalyzer(nm_var.split(":")[0], nm_var.split(".")[-1], "c."):
						warnings.warn("Mutalyzer failed at " + str(nm_var) + "\t" + str(nc_var))
						error = ("Mutalyzer failed  at " + str(nm_var) + "\t" + str(nc_var))
					else:
						warnings.warn("Mismatch at " + str(nm_var) + "\t" + str(nc_var))
						error = ("Mismatch at " + str(nm_var) + "\t" + str(nc_var))
				nc_pos = tools.get_num(nc_var.split(":")[-1])
				realpos = nc_pos-minpos+250
				wtstr = seq[realpos - 250:realpos + 251]
				mutstr = tools.update_str(seq, nc_change[-1], realpos)[realpos - 250:realpos + 251]
				if realpos > len(seq) or seq[realpos] != nc_change[0]:
					warnings.warn("NTs don't match at " + str(nm_var) + "\t" + str(nc_var) + "\n" + str(wtstr))
					error = ("NTs don't match at " + str(nm_var) + "\t" + str(nc_var) + "\n" + str(wtstr))
				else:
					error = ""
				if minus:
					wtstr = "".join([complement[wtstr[i]] for i in range(len(wtstr))])
					mutstr = "".join([complement[mutstr[i]] for i in range(len(mutstr))])
				if flip:
					wtstr=wtstr[::-1]
					mutstr=mutstr[::-1]
				muts[nm_var] = {"WT":wtstr, "MUT":mutstr, "startpos":nc_pos - 250, "endpos":nc_pos + 250, "error":error, "mixup":mixup, "minus":minus, "flip":flip}
				if error != "":
					try:
						print(nm_var)
						print(str(z) + "\t" + nc_var + "\n" + wtstr[250-20:250] + " " + wtstr[250] + " " + wtstr[251:250+20])
					except IndexError as e:
						print(e)
			except Exception as e:
				print(str(nm_var) + "\t" + str(e))
	return(muts)

def query_gscholar(var):
	i = 0
	html = []
	while True:
		try:
			r = tools.retry_request(requests.get, ["https://scholar.google.com/scholar?start=" + str(i*10) + "&q=" + urllib.parse.quote(var)])
			pieces = re.split("\<div class\=\"gs_r gs_or gs_scl\"", r.text)
			if len(pieces) > 1:
				pieces = pieces[1:]
			else:
				break
			for piece in pieces:
				links = re.findall("href\=\"(https?\:\/\/.+?)\"", piece)
				if len(links) > 0:
					html.append(links[0])
			i = i + 1
		except Exception as e:
			break
	return(html)

def get_all_variant_info(ncvar):
	data = query_clinvar(ncvar)
	testlist = [ncvar]
	assembly = [k for k in chromosome_map.keys() if ncvar.split(":")[0] in chromosome_map[k].values()][0]
	if "rsid" in data:
		freqs, clinsig, nm, rsid = query_dbSNP(data["rsid"],assembly)
		if not data["rsid"].startswith("rs"):
			testlist.append("rs" + data["rsid"])
		else:
			testlist.append(data["rsid"])
	else:
		freqs, clinsig, nm, rsid = query_dbSNP(ncvar,assembly)
		data["rsid"] = rsid
		testlist.append(rsid)
	data.update(freqs)
	if "NM_acc" in data:
		testlist += data["NM_acc"]
	data["google"] = []
	for acc in testlist:
		data["google"] += query_gscholar(acc)
	data["google"] = list(set(data["google"]))
	if "disease" in data:
		data["disease"] += clinsig
	else:
		data["disease"] = clinsig
	if "NM_acc" in data:
		data["NM_acc"] += nm
	else:
		data["NM_acc"] = nm
	data["NM_acc"] = [max(re.findall("NM\_\d+\.?\d*.*?\:c\.[\-|\+|\*]?\d*[\-|\+|\*]?\d+[A-Z]+\>[A-Z]+", nmacc), key=len) for nmacc in data["NM_acc"]]
	return(data)

def query_clinvar(ncvar):
	data = {}
	r = tools.retry_func(Entrez.esearch, [], {'db':'clinvar', 'term':ncvar}).read()
	ids = re.findall("\<Id\>(\d+)\<\/Id\>", r)
	for idi in ids:
		#handle = tools.retry_func(Entrez.efetch, [], {'db':"clinvar", 'is variation': True, 'id':idi, 'rettype':"vcv", 'retmode':"text", "from esearch":True})
		r = tools.retry_request(requests.get, ['https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id=' + idi + '&from_esearch=true'], {})
		disease = re.findall("\<RCVAccession Title\=.+Interpretation=\"(.+?)\".+\>", r.text)
		if disease != None and len(disease) > 0:
			data['disease'] = disease
		rsid = list(set(re.findall("\<XRef .*?ID\=\"\d+\".*?\>", r.text)))
		rsid = [max(re.findall("\<XRef .*?ID\=\"(\d+)\".*?\>", rsid[i]), key=len) for i in range(len(rsid)) if "dbSNP" in rsid[i] or "rs" in rsid[i]]
		if rsid != None and len(rsid) > 0:
			data['rsid'] = max(set(rsid), key=rsid.count)
		nm_id = list(set(re.findall("\<Name\>(NM\_\d+\.?\d*.*?\:c\.[\-|\+|\*]?\d*[\-|\+|\*]?\d+[A-Z]+\&gt\;[A-Z]+).*\<\/Name\>", r.text)))
		nm_id = [nm_id[i].replace("&gt;", ">") for i in range(len(nm_id))]
		if len(nm_id) > 0:
			gene_names = [re.findall("\(([A-Za-z\d]+)\)", nm_id[i]) for i in range(len(nm_id))]
			gene_names = [i for j in gene_names for i in j]
			nm_id = [re.sub("\([A-Za-z\d]+\)", "", nm_id[i]) for i in range(len(nm_id))]
			data["gene name"] = gene_names
			data["NM_acc"] = nm_id
	return(data)

def query_dbSNP(rsid, altonly=True, refineorder=['gnomAD - Genomes', '1000Genomes', 'HapMap', 'gnomAD - Exomes'], verbose=False):
	if not (tools.is_numeric(rsid) or (isinstance(rsid, str) and rsid.startswith("rs"))):
		if "nc" in rsid.lower():
			try:
				assembly = [k for k in chromosome_map.keys() if rsid.split(":")[0] in chromosome_map[k].values()][0]
				ch = [k for k in chromosome_map[assembly].keys() if chromosome_map[assembly][k] == rsid.split(":")[0]][0]
				base_pos = tools.get_num(rsid.split(".")[-1])
				if assembly == "GRCh37":
					r = tools.retry_request(requests.get, ["https://www.ncbi.nlm.nih.gov/snp/?term=(" + str(ch) + "%5BChromosome%5D)%20AND%20" + str(base_pos) + "%5BBase%20Position%20Previous%5D"])
				else:
					r = tools.retry_request(requests.get, ["https://www.ncbi.nlm.nih.gov/snp/?term=(" + str(ch) + "%5BChromosome%5D)%20AND%20" + str(base_pos) + "%5BBase%20Position%5D"])
				ids = re.findall("rs\d+", r.text)
				rsid = ids[0]
			except Exception as e:
				return({})
	rsid = rsid.replace("rs", "")
	r = tools.retry_request(requests.get, ["https://www.ncbi.nlm.nih.gov/snp/rs" + str(rsid)], {})
	nm = re.findall("NM\_\d+\.?\d*.*?\:c\.[\-|\+|\*]?\d*[\-|\+|\*]?\d+[A-Z]+\&gt\;[A-Z]+", r.text) + re.findall("NM\_\d+\.?\d*.*?\:c\.[\-|\+|\*]?\d*[\-|\+|\*]?\d+[A-Z]+\>[A-Z]+", r.text)
	nm = list(set([nmi.replace("&gt;", ">") for nmi in nm]))
	parsed_html = BeautifulSoup(r.text, features="lxml")
	freqs = parsed_html.body.find_all("tr", attrs={"class":"chi_row"}) + parsed_html.body.find_all("tr", attrs={"class":"par_row"})
	freqd = {}
	alt = ""
	for i, freq in enumerate(freqs):
		names = tuple(re.findall("\<a.+\>(.+?)\<\/a\>", str(freq)))
		alleles = re.findall("\<td\>([A-Za-z]+\=\d?\.?\d*)\<\/td\>", str(freq))
		alleles = {max(re.findall("[A-Za-z]+", allele), key=len):tools.get_num(allele) for allele in alleles}
		if altonly: #Only give frequency for ALT allele
			if len(alleles.keys()) < 2:
				continue
			elif alt == "" and len(alleles.keys()) > 1:
				alt = list(alleles.keys())[1]
			elif alt != list(alleles.keys())[1]:
				warnings.warn("Not all alleles are ordered the same")
			alleles = alleles[alt]
		freqd[names] = alleles
	for order in refineorder: #Include results for only one study, depending on which studies have data
		if any([True for k in freqd if order == k[0]]):
			freqd = {k[1]:v for k,v in freqd.items() if order == k[0]}
			if verbose:
				print("Using study " + order)
			break;
	clinsig = []
	clin = parsed_html.body.find_all("div", attrs={"id":"clinical_significance"})
	if len(clin) > 0:
		clin = clin[0]
		alleles = re.split("\<div class\=\"sect\_heading\"\>", str(clin))
		for allele in alleles:
			nt = re.findall("Allele:\s*([A-Z])", allele)
			if len(nt) < 1 or nt[0] != alt:
				continue
			nt = nt[0]
			rows = allele.split("<tr")
			for row in rows:
				cols = row.split("<td>")
				cols = [col.split("</td>")[0] for col in cols]
				if len(cols) >= 4:
					clinsig.append((cols[2], cols[3]))
	return(freqd, clinsig, nm, rsid)

def compare_splicing(wtstr, mutstr, loc=None, genename="pholder", directory="./sources/", normalize = False):
	if not os.path.isdir(directory):
		os.mkdir(directory)
	data = {}
	if loc == None:
		loc = len(wtstr)//2
		if not len(wtstr) % 2 == 1:
			warnings.warn("Location not given")
	wt_esefinder = compute_features.run_ESEfinder_wrapper(genomic_aligned=wtstr, ORF_aligned=wtstr, genename=genename, space="-", maxlen=4900, quiet=False)
	mut_esefinder = compute_features.run_ESEfinder_wrapper(genomic_aligned=mutstr, ORF_aligned=mutstr, genename=genename, space="-", maxlen=4900, quiet=False)
	data["esefinder"] = any([wt_esefinder[i] != mut_esefinder[i] for i in range(len(wt_esefinder))])
	data["fas_ess"] = compute_features.fas_ess(wtstr) == compute_features.fas_ess(mutstr)
	data["exonscan"] = compute_features.get_exonscan(wtstr) == compute_features.get_exonscan(mutstr)
	for f in ["ESR.tsv", "Z_EI.tsv", "Z_WS.tsv"]:
		df = pd.read_csv(os.path.join(directory, f), sep="\t", index_col=0, header=None).T.to_dict()
		df = {k:v[1] for k,v in df.items()}
		if normalize:
			u = statistics.mean(df.values())
			o = statistics.stdev(df.values())
			#df = {k:(v-u)/o for k,v in df.items()}
		wtscore = []
		mutscore = []
		for i in range(6):
			wtscore.append(df[wtstr[loc-i:loc-i+6]])
			mutscore.append(df[mutstr[loc-i:loc-i+6]])
		wtscore = sum(wtscore)
		mutscore = sum(mutscore)
		data[f.split(".")[0]] = mutscore - wtscore
	return(data)

def seq_loglikelihood(alignment, focus=None, positions=None, space="-", verbose=False):
	if focus == None or focus not in alignment:
		focus = list(alignment.keys())[0]
	if positions == None:
		positions = np.arange(0, len(alignment[focus].replace(space, "")))
	dist = score_variants.compute_distribution(alignment, focus)
	converted_pos = []
	for pos in positions:
		i = compute_features.convert_position2(alignment[focus], pos)
		converted_pos.append(i)
	loglikelihood = {}
	matching = {}
	tmpseq = "".join([alignment[focus][pos] for pos in converted_pos])
	for k,v in alignment.items():
		loglikelihood[k] = 0
		matching[k] = 0
		good_positions = 0
		for i, pos in enumerate(converted_pos):
			try:
				loglikelihood[k] += math.log(dist[i+1][v[pos]])
				if v[pos] == alignment[focus][pos]:
					matching[k] += 1
				elif k == focus:
					print(v[pos] + "\t" + alignment[focus][pos])
				good_positions += 1
			except (KeyError, ValueError) as e:
				if verbose:
					print(str(k) + "\t" + str(i) + "\t" + str(e))
		if good_positions != 0:
			loglikelihood[k] /= good_positions
			matching[k] /= len(converted_pos)
		elif loglikelihood[k] != 0:
			warnings.warn("Loglikelihood is nonzero, but number of positions is not")
	normal = loglikelihood[focus]
	loglikelihood = {k:v-normal for k,v in loglikelihood.items()}
	return(loglikelihood, matching)

def parse_miRDB_info(f, regex="\d.NC\_\d+\.\d+.g\.\d+[ACGT].[ACGT]", maxval = 100):
	data = {}
	wb = xlrd.open_workbook(f)
	for i, sheetname in enumerate(wb.sheet_names()):
		#print(f + "\t" + sheetname)
		sheetname = str(sheetname)
		sheet = wb.sheet_by_index(i)
		loc = {"WT":None, "MUT":None}
		wtloc = None
		mutloc = None
		if re.search(regex, sheetname):
			data[sheetname] = {"WT":{}, "MUT":{}}
			for k in range(0, maxval):
				for l in range(0, maxval):
					try:
						if str(sheet.cell_value(k, l)).strip() == "MUT" and loc["MUT"] == None:
							loc["MUT"] = (k+2, l)
							for m in range(max(0, l-3), min(maxval, l+3)):
								if str(sheet.cell_value(k+1, m)).strip() == "Target Rank":
									loc["MUT"] = (k+2, m)
									break
						elif str(sheet.cell_value(k, l)).strip() == "WT" and loc["WT"] == None:
							loc["WT"] = (k+2, l)
							for m in range(max(0, l-3), min(maxval, l+3)):
								if str(sheet.cell_value(k+1, m)).strip() == "Target Rank":
									loc["WT"] = (k+2, m)
									break
						#print(str((k,l)) + "\t" + str(sheet.cell_value(k,l)))
					except (IndexError, ValueError) as e:
						pass
			if loc["WT"] == None or loc["MUT"] == None:
				print(sheetname + "\t" + "MUT and WT row not found: " + str(sheetname) + "\t" + str(loc))
			#elif loc["WT"][0] != loc["MUT"][0]:
			#	print(sheetname + "\t" + "MUT and WT rows not same: " + str(loc))
			else:
				for k in range(loc["MUT"][0], loc["MUT"][0] + maxval):
					for key in loc.keys():
						try:
							rank = int(sheet.cell_value(k, loc[key][1]))
							score = int(sheet.cell_value(k, loc[key][1] + 1))
							miRNA = str(sheet.cell_value(k, loc[key][1] + 2))
							data[sheetname][key][miRNA] = {"Rank":rank, "Score":score}
						except (IndexError, ValueError) as e:
							pass
				for k in data[sheetname].keys():
					if data[sheetname][key]  == {}:
						print(f + "\t" + sheetname + "\t" + key + "\t" + str(loc[key]))
	return(data)

#finds the sequences and accession ids for a given gene
#completely finds and aligns sequences for a gene in preparation for compute_features.precompute
def pipeline(genename, outdir="./genes/", organism="homo sapiens", write=False, email=Entrez.email):
	if not os.path.isdir(outdir):
		os.mkdir(outdir)
	Entrez.email=email
	seqids = get_accids(genename, organism=organism)
	seqs, cdses, mRNA, protein = get_seqs_and_exons(genename, seqids["genomic"], seqids["mRNA"], organism=organism)
	if write:
		tools.write_fasta(seqs, os.path.join(outdir, genename+".fasta"))
		tools.write_fasta(cdses, os.path.join(outdir, genename+".cds"))
		parse_range.parse_range(os.path.join(outdir, genename + ".cds"))
	else:
		return(seqs, cdses, mRNA, protein)
