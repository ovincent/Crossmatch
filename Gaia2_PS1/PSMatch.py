import photConv as photConv
from astroquery.vizier import Vizier #bugged in astropy 4.0+
import astropy.units as unitz
import astropy.coordinates as coord

v = Vizier(columns=['DR2Name', 'RA_ICRS','DE_ICRS', 'Gmag', 'BP-RP', 'Plx', 'e_Plx'])

def find_ps1(gaiainfo=[],SDSSPhot=[],searchR=0.00139):
    uSDSS,gSDSS,rSDSS,iSDSS,zSDSS=SDSSPhot
    gaianame,RA,DEC,G,BR,p=gaiainfo
    
    flags,info=[False,'-'],['-']*13
    
    #convert epoch here
    #...
    
    #cherche par coord avec searchR sec d'arc
    results=v.query_region(coord.SkyCoord(ra=RA, dec=DEC,unit=(unitz.hour, unitz.deg),frame='icrs'),
                         radius=coord.Angle(searchR,'deg'),
                         catalog='II/349/ps1')
    if len(results)>0: flags,info=PS1search(results,RA,DEC,G,BR,gSDSS)
    
    #plus grand rayon
    else:
        results=v.query_region(coord.SkyCoord(ra=RA, dec=DEC,unit=(unitz.hour, unitz.deg),frame='icrs'),
                         radius=coord.Angle(searchR*3,'deg'),
                         catalog='II/349/ps1')
        flags,info=PS1search(results,RA,DEC,G,BR,gSDSS)
    
    return flags,info

def PS1search(table,RA,DEC,G,BR,gSDSS,maxMagDiff=0.3):
    table=table[0]
    
    #default values
    ps1matched=False
    ps1meth='-'
    name,ra,dec,g,ge,r,re,i,ie,z,ze,y,ye=['-']*13
    
    RA_all=table['RAJ2000']
    DEC_all=table['DEJ2000']
    g_all=table['gmag']
    
    #si on a un g
    if gSDSS!='00.00':
        #cherche l'objet avec le g le plus proche
        idx=(np.abs(g_all - gSDSS)).argmin()
        if np.abs(gSDSS-g_all[idx])<maxMagDiff:
            ps1matched=True
            ps1meth='gSDSS'
            name,ra,dec=table['objID'][idx],table['RAJ2000'][idx],table['DEJ2000'][idx]
            g,ge=table['gmag'][idx],table['e_gmag'][idx]
            r,re=table['rmag'][idx],table['e_rmag'][idx]
            i,ie=table['imag'][idx],table['e_imag'][idx]
            z,ze=table['zmag'][idx],table['e_zmag'][idx]
            y,ye=table['ymag'][idx],table['e_ymag'][idx]
    
    #sinon on convertit G gaia en g SDSS
    if not ps1matched:
        g_approx=photConv.gaia2SDSS(G,BR)
        idx=(np.abs(g_all - g_approx)).argmin()
        if np.abs(g_approx-g_all[idx])<maxMagDiff:
            ps1matched=True
            ps1meth='GaiaG'
            name,ra,dec=table['objID'][idx],table['RAJ2000'][idx],table['DEJ2000'][idx]
            g,ge=table['gmag'][idx],table['e_gmag'][idx]
            r,re=table['rmag'][idx],table['e_rmag'][idx]
            i,ie=table['imag'][idx],table['e_imag'][idx]
            z,ze=table['zmag'][idx],table['e_zmag'][idx]
            y,ye=table['ymag'][idx],table['e_ymag'][idx]
    
    return [ps1matched,ps1meth],[name,ra,dec,g,ge,r,re,i,ie,z,ze,y,ye]