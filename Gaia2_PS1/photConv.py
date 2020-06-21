def BVRI2SDSS(BVRI=[]):
	"""
	Johnson to SDSS. Converts B to g.
	based on https://arxiv.org/pdf/astro-ph/0609121.pdf
	Input: BVRI photometry array or list (in this order).
	Return: g photometry in SDSS system.
	"""
	return BVRI[0] -0.370*(BVRI[0]-BVRI[1]) -0.124

def BVRI2gaia(BVRI=[]):
	"""
	BVRI to Gaia. Converts B-V to G.
	Valid for -0.3 < B-V < 2.4
	based on https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_PhotTransf.html
	Input: BVRI photometry array or list (in this order).
	Return: G photometry in Gaia system.
	"""
	return BVRI[1] -0.02907 -0.02385*(BVRI[0]-BVRI[1]) -0.2297*(BVRI[0]-BVRI[1])**2 -0.001768*(BVRI[0]-BVRI[1])**3

def JHK2gaia(JHK=[]):
	"""
	Johnson to Gaia. Converts H-K to G.
	Valid for −0.3 < H−K < 2.2
	based on https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_PhotTransf.html
	Input: JHK photometry array or list (in this order).
	Return: G photometry in Gaia system.
	"""
	return JHK[2]+ 0.23587+ 4.0548*(JHK[0]-JHK[2]) -2.5608*(JHK[0]-JHK[2])**2 +2.2228*(JHK[0]-JHK[2])**3 -0.54944*(JHK[0]-JHK[2])**4
	
def gaia2SDSS(G,BR):
	"""
	Gaia to SDSS. Converts G and B-R to g.
	Valid for −0.3 < H−K < 2.2
	based on https://arxiv.org/pdf/1804.09368.pdf
	Input: G and B-R Gaia photometry.
	Return: g photometry in SDSS system.
	"""
	return G -0.13518 +0.46245*BR +0.25171*BR**2 -0.021349*BR**3