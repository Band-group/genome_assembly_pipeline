def find_input_fastqs( reads, sample, type, prep = None, replicate = None, miniaturisation = None ):
	if type == 'pacbio':
		if prep is None:
			return reads.loc[ reads['sample'] == sample ].filename.tolist()
		else:
			result = (
				reads.loc[
					(reads['sample'] == sample)
					& (reads['prep'] == prep)
					& (reads['replicate'] == int(replicate))
					& (reads['miniaturisation'] == miniaturisation)
				].filename.tolist()
			)
			return result
	elif type == 'illumina':
		return config['illumina-reads'][sample]
	elif type == 'both':
		return (config['illumina-reads'][sample] + reads.loc[ reads['sample'] == sample ].filename.tolist())