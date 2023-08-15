#########################################################################################################
import pandas as pd
import matplotlib.pyplot as plt
#import netCDF4 as nc
import numpy as np
import geopandas as gpd
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union, cascaded_union

from matplotlib import colors
import seaborn                   # Graphics
import geopandas as gpd          #street_typ Spatial data manipulation
import pandas                    # Tabular data manipulation
import rioxarray                 # Surface data manipulation
import xarray                    # Surface data manipulation

from pysal.explore import esda   # Exploratory Spatial analytics
from pysal.lib import weights    # Spatial weights
#import contextily                # Background tiles

import tkinter
import matplotlib
matplotlib.use('TkAgg')

#########################################################################################################

def calc_gord(g,db):
    # Break observations into significant or not
    sig = g.p_sim < 0.025
    # Plot non-significant clusters
    ns = db.loc[sig==False, 'geometry']
    nsu = db.loc[sig==False, 'u']
    # Plot HH clusters
    hh = db.loc[(g.Zs > 0) & (sig==True), 'geometry']
    hu = db.loc[(g.Zs > 0) & (sig==True), 'u']
    # Plot LL clusters
    ll = db.loc[(g.Zs < 0) & (sig==True), 'geometry']
    lu = db.loc[(g.Zs < 0) & (sig==True), 'u']
    # set up db
    classes = ['high']*len(hh)+['low']*len(ll)+['ns']*len(ns)
    geo = list(hh)+list(ll)+list(ns)
    d = pd.DataFrame(classes)
    d['geometry'] = geo
    d['u'] = list(hu)+list(lu)+list(nsu)
    gd = gpd.GeoDataFrame(d,crs = db.crs)
    gd.columns = ['getisord','geometry','u']
    return gd


def get_ordg_hotspots_from_shp(gdf,w_surface,v):
    ordg = []
    for i in range(len(v)):
        # Gi
        go_i = esda.getisord.G_Local(gdf[v[i]], w_surface)
        ordg.append(calc_gord(go_i,gdf))
    #
    # clean up and reocmbine
    for i in range(len(v)):
        gdf = pd.merge(ordg[i],gdf,on='u')
        gdf['getisord_'+v[i]] = gdf.getisord
        gdf['geometry'] = gdf.geometry_x
        gdf = gdf.drop('geometry_y',axis=1)
        gdf = gdf.drop('geometry_x',axis=1)
        gdf = gdf.drop('getisord',axis=1)
    return gdf


def get_ordg_hotspots_from_shp2(gdf,w_surface,v):
    ordg = []
    for i in range(len(v)):
        # Gi
        go_i = esda.getisord.G_Local(gdf[v[i]], w_surface)
        # Break observations into significant or not
        #
        g,db = go_i,gdf
        sig = g.p_sim < 0.05/3
        # Plot non-significant clusters
        ns = db.loc[sig==False, ['geometry','u']]
        ns['getisord_'+v[i]] = 'ns'
        #nsu = db.loc[sig==False, 'u']
        # Plot HH clusters
        hh = db.loc[(g.Zs > 0) & (sig==True), ['geometry','u']]
        hh['getisord_'+v[i]] = 'high'
        #hu = db.loc[(g.Zs > 0) & (sig==True), 'u']
        # Plot LL clusters
        ll = db.loc[(g.Zs < 0) & (sig==True), ['geometry','u']]
        ll['getisord_'+v[i]] = 'low'
        #lu = db.loc[(g.Zs < 0) & (sig==True), 'u']
        # set up db
        d = hh.append(ns).append(ll)
        gd = gpd.GeoDataFrame(d,crs = db.crs)
        ordg.append(gd)
    #
    # clean up and reocmbine
    for i in range(len(v)):
        gdf = pd.merge(ordg[i],gdf,on='u')
        #gdf['getisord_'+v[i]] = gdf.getisord
        gdf['geometry'] = gdf.geometry_x
        gdf = gdf.drop('geometry_y',axis=1)
        gdf = gdf.drop('geometry_x',axis=1)
        #gdf = gdf.drop('getisord',axis=1)
    return gdf

#########################################################################################################
# Pull in files
# for plotting
streets = gpd.read_file('/projects/b1045/montgomery/shapefiles/geo_export_28abe032-5cb4-4d93-9c9e-cf51d26e8169.shp')
class1 = streets[streets['class'] == '1'].reset_index(drop=True)
class2 = streets[streets['class'] == '2'].reset_index(drop=True)


#
demos = gpd.read_file('/projects/b1045/montgomery/Paper2/demo_chi_of_baseline.shp') # double counts hisplatinx
demos = demos.drop(['White','Black', 'Latino', 'NativeAmer', 'Asian', 'HawaainPac', 'Other', 'Two+',
    'pWhite', 'pBlack','pLatino', 'pNativeAme', 'pAsian', 'pHawaainPa', 'pOther', 'pTwo+','GISJOIN','geometry'],axis=1)

shp_tracts = gpd.read_file('/projects/b1045/montgomery/Paper2/demo_chi_of_baseline_hisplatinx.shp')
shp_tracts.columns = ['total', 'total_NHL', 'White', 'Black', 'Native', 'Asian',
       'HawaiiPI', 'Other', 'Two', 'ALUKE010', 'ALUKE011', 'HispLatino',
       'GISJOIN', 'u', 'geometry']

demos = pd.merge(shp_tracts,demos,on=['u'])
demos['pWhite'] = demos.White/demos.total
demos['pBlack'] = demos.Black/demos.total
demos['pHispLatino'] = demos.HispLatino/demos.total
demos['pAsian'] = demos.Asian/demos.total

# file just made
dg = gpd.read_file('/projects/b1045/montgomery/microsoft/all_data_grid_eclipse_rf.geojson')

# need to crop final datasets to chicago limits
from shapely.ops import cascaded_union, unary_union
chi_shp = gpd.read_file('/projects/b1045/montgomery/shapefiles/Chicago/commareas/geo_export_77af1a6a-f8ec-47f4-977c-40956cd94f97.shp')
chi_shp = chi_shp[chi_shp['community'] != "OHARE"].reset_index()
chi_outer = cascaded_union(chi_shp.geometry)
# clip to chicago only
dg = gpd.clip(dg,chi_outer)
demos = demos.to_crs(dg.crs)
#demos = gpd.clip(demos,chi_outer)
gp_chi_outer = gpd.GeoDataFrame({'geometry':chi_outer},crs=chi_shp.crs)

#
ndvifeb = np.array(pd.read_csv('NDVI_crop_feb2022.csv',index_col=0)).ravel()
ndviaug = np.array(pd.read_csv('NDVI_crop_aug2021.csv',index_col=0)).ravel()

fin = pd.read_csv('reanalysis_feb2022.csv',index_col = 0)
fin.columns = ['u','winddir_rean02','airt_rean02','sf_rean02','rh_rean02','ws_rean02','hr']

d ='/projects/b1045/montgomery/microsoft/'
tmp_shp = gpd.read_file('/projects/b1045/montgomery/microsoft/rectangle_chicago_1km_grid.shp')
tmp_shp['ndvi_aug'] = ndviaug
tmp_shp['ndvi_feb'] = ndvifeb
tmp_shp['trop_02'] = np.array(pd.read_csv(d+'02_2022_tropomi_NO2.csv',index_col = 0)).ravel()
tmp_shp['trop_08'] = np.array(pd.read_csv(d+'08_2021_tropomi_NO2.csv',index_col = 0)).ravel()
tmp_shp['trop_02_ct'] = np.array(pd.read_csv(d+'pixel_counts_202108.csv',index_col = 0)).ravel()
tmp_shp['trop_08_ct'] = np.array(pd.read_csv(d+'pixel_counts_202202.csv',index_col = 0)).ravel()
tmp_shp['cmaq_08'] = np.array(pd.read_csv(d+'cropped_no2_monthly_ESPG4326.csv',index_col = 0)).ravel()
tmp_shp['cmaq_02'] = np.array(pd.read_csv(d+'cropped_no2_monthly_feb2022_ESPG4326.csv',index_col = 0)).ravel()
tmp_shp['trop_count_08']= np.array(pd.read_csv(d+'pixel_counts_202108.csv',index_col = 0)).ravel()
tmp_shp['trop_count_02']= np.array(pd.read_csv(d+'pixel_counts_202202.csv',index_col = 0)).ravel()
tmp_shp = tmp_shp.drop('geometry',axis=1)
tmp_shp = pd.merge(tmp_shp,fin,on='u')

dg = gpd.GeoDataFrame(pd.merge(tmp_shp,dg,on='u'),crs = dg.crs)

#dg['geometry'] = dg.geometry_y
#dg=dg.drop('geometry_y',axis=1)
#dg=dg.drop('geometry_x',axis=1)

dg = pd.merge(dg,demos,on='u')
dg['geometry'] = dg.geometry_x
dg=dg.drop('geometry_y',axis=1)
dg=dg.drop('geometry_x',axis=1)
#dg = gpd.clip(dg,chi_outer)

v = ['cmaq_08','trop_08','rf_month_august']
v2 = ['cmaq_02','trop_02','rf_month_february']

#make arrays
w_surface = weights.distance.Kernel.from_dataframe(dg, fixed=False, k=4)

w_surface = weights.Queen.from_dataframe(dg)
dfa = get_ordg_hotspots_from_shp2(dg,w_surface,v)
dgs = gpd.GeoDataFrame(dfa,crs="EPSG:4326")

# Same for febrya
dgwa = get_ordg_hotspots_from_shp2(dg,w_surface,v2)
dgw = gpd.GeoDataFrame(dgwa,crs="EPSG:4326")

#########################################################################################################
#plot raw data normalized
import string
alphabet = np.array([string.ascii_lowercase[i] for i in range(26)])

def norm_data(x):
    return (x - np.quantile(x,0.01))/(np.quantile(x,0.99)-np.quantile(x,0.01))

dgs['ec_08_norm'] = norm_data(dgs.rf_month_august)
dgs['cmaq_08_norm'] = norm_data(dgs.cmaq_08)
dgs['trop_08_norm'] = norm_data(dgs.trop_08)
dgw['ec_02_norm'] = norm_data(dgw.rf_month_february)
dgw['cmaq_02_norm'] =  norm_data(dgw.cmaq_02)
dgw['trop_02_norm'] = norm_data(dgw.trop_02)

# plot each cluster
fig,axs = plt.subplots(2,3,figsize=(8,6))
axs=axs.ravel()
titles=['Eclipse','WRF-CMAQ','TropOMI','','','']
labels = ['August 2021','','','February 2022','','']
di = ['ec_08_norm','cmaq_08_norm','trop_08_norm','ec_02_norm','cmaq_02_norm','trop_02_norm']

for i in range(len(axs)):
    ax=axs[i]
    dgw.plot(facecolor='lightgray',ax=ax)
    if i < 3:
        dgs.plot(di[i],ax=ax,vmin = 0,vmax = 1,legend=False)
        ax.text(0.25, 0.03, '$\mu$ = %.1f'%(dgs[di[i]].mean()), va='bottom', ha='center',transform=ax.transAxes,fontsize = 9,c='k')
    else: 
        dgw.plot(di[i],ax=ax,vmin = 0,vmax = 1,legend=False)
        ax.text(0.25, 0.03, '$\mu$ = %.1f'%(dgw[di[i]].mean()), va='bottom', ha='center',transform=ax.transAxes,fontsize = 9,c='k')
    #if i == 4: dgw.plot(di[i],ax=ax,vmin = 0,vmax = 1,legend=True,legend_kwds={"shrink":.5,'orientation': "horizontal"})
    ax.axis('off')
    ax.set_title(titles[i])
    ax.text(-0.03, 0.9, '(%s)'%(alphabet[i]), va='bottom', ha='center',transform=ax.transAxes,fontsize = 9,c='k')
    ax.text(-0.05, 0.3, '%s'%(labels[i]), va='bottom', ha='center',transform=ax.transAxes,fontsize = 9,c='k',rotation=90)

norm = Normalize(vmin=0, vmax=1)
n_cmap = cm.ScalarMappable(norm=norm, cmap="viridis")
n_cmap.set_array([])
cbar = ax.get_figure().colorbar(n_cmap, ax=axs, orientation='vertical',shrink = 0.5,label=r'Normalized NO$_2$')
cbar.ax.tick_params(labelsize=9)

plt.savefig('normalized_datasets.png',dpi=300)
plt.show()

# plot each cluster
fig,axs = plt.subplots(2,3)
axs=axs.ravel()
titles=['(a) G* CMAQ','(b) G* TropOMI','(c) G* Eclipse','','','']
datas = [dgs,dgs,dgs,dgw,dgw,dgw]
vs = v+v2
axs = axs.ravel()
for i in range(len(vs)):
    ax=axs[i]
    dgw.plot(facecolor='lightgray',ax=ax)
    datas[i][datas[i]['getisord_'+vs[i]] == 'high'].plot(ax=ax,facecolor='red',alpha=0.8,label='High August',legend=True)
    #dgw[dgw['getisord_'+v[i]] == 'low'].plot(ax=ax,facecolor='mediumblue',alpha=0.8,label='Low August',legend=True)
    plt.rcParams['hatch.linewidth'] = 0.5
    ax.set_title(titles[i])
   # dgw[dgw['getisord_'+vs[i]] == 'high'].plot(ax=ax,facecolor='red',alpha=0.3,hatch='//',label='High February',legend=True)
    #dgw[dgw['getisord_'+v2[i]] == 'low'].plot(ax=ax,facecolor='mediumblue',alpha=0.3,hatch='//',label='Low February',legend=True)
    class1.plot(facecolor='None',edgecolor='k',alpha=0.8,linewidth=0.8,ax=ax)
    ax.axis('off')


plt.tight_layout()
plt.savefig('/projects/b1045/montgomery/Paper2/clustering_on_each_dataset_2row.png',dpi=300)
plt.show()



#########################################################################################################
# Seperate Groups and add labels
dgw['HH'] = '' #eclipse + tropomi + cmaq
dgs['LL'] = ''
dgs['SW'] = '' # eclipse + tropomi
dgs['CL2'] = '' # eclipse + cmaq
dgs['TC'] = '' # tropomi + cmaq
dgs['CL'] = '' #cmaq only
dgs['EO'] = '' #eclipse only
dgs['TO'] = '' #tropomi only

dgs.loc[(dgs.getisord_cmaq_08 == 'high') & (dgs.getisord_trop_08 == 'high') & (dgs.getisord_rf_month_august == 'high'),'HH'] = 'HH'
dgs.loc[(dgs.getisord_cmaq_08 == 'low') & (dgs.getisord_trop_08 == 'low') & (dgs.getisord_rf_month_august == 'low'),'LL'] = 'LL'

# winter
dgw['HH'] = ''
dgw['LL'] = ''
dgw['SW'] = ''
dgw['CL'] = ''
dgw['CL2'] = ''
dgw['EO'] = '' #eclipse only
dgw['TO'] = '' #tropomi only
dgw['TC'] = '' # tropomi + cmaq

dgw.loc[(dgw.getisord_cmaq_02 == 'high') & (dgw.getisord_trop_02 == 'high') & (dgw.getisord_rf_month_february == 'high'),'HH'] = 'HH'
dgw.loc[(dgw.getisord_cmaq_02 == 'low') & (dgw.getisord_trop_02 == 'low') & (dgw.getisord_rf_month_february == 'low'),'LL'] = 'LL'

#########################################################################################################
# LOW CONFIDENCE

dgs.loc[((dgs.getisord_cmaq_08 == 'high') & 
    ((dgs.getisord_trop_08 == 'low') | (dgs.getisord_trop_08 == 'ns')) 
    & ((dgs.getisord_rf_month_august == 'low') | (dgs.getisord_rf_month_august == 'ns'))),'CL'] = 'CL'


dgs.loc[((dgs.getisord_cmaq_08 == 'high') & 
    ((dgs.getisord_trop_08 == 'low') | (dgs.getisord_trop_08 == 'ns')) 
    & ((dgs.getisord_rf_month_august == 'high') | (dgs.getisord_rf_month_august == 'high'))),'CL2'] = 'CL2'


dgw.loc[((dgw.getisord_cmaq_02 == 'high') & 
    ((dgw.getisord_trop_02 == 'low') | (dgw.getisord_trop_02 == 'ns')) 
    & ((dgw.getisord_rf_month_february == 'low') | (dgw.getisord_rf_month_february == 'ns'))),'CL'] = 'CL'

dgw.loc[((dgw.getisord_cmaq_02 == 'high') & 
    ((dgw.getisord_trop_02 == 'low') | (dgw.getisord_trop_02 == 'ns')) 
    & ((dgw.getisord_rf_month_february == 'high') | (dgw.getisord_rf_month_february == 'high'))),'CL2'] = 'CL2'

dgs.loc[((dgs.getisord_trop_08 == 'high') & (dgs.getisord_rf_month_august == 'high')) & 
    (((dgs.getisord_cmaq_08 == 'low') | (dgs.getisord_cmaq_08 == 'ns'))),'SW'] = 'SW'

dgw.loc[((dgw.getisord_trop_02 == 'high') & (dgw.getisord_rf_month_february == 'high')) & 
    (((dgw.getisord_cmaq_02 == 'low') | (dgw.getisord_cmaq_02 == 'ns'))),'SW'] = 'SW'


dgs.loc[((dgs.getisord_cmaq_08 == 'high') & 
    ((dgs.getisord_trop_08 == 'high') | (dgs.getisord_trop_08 == 'high')) 
    & ((dgs.getisord_rf_month_august == 'low') | (dgs.getisord_rf_month_august == 'ns'))),'TC'] = 'TC'

dgw.loc[((dgw.getisord_cmaq_02 == 'high') & 
    ((dgw.getisord_trop_02 == 'high') | (dgw.getisord_trop_02 == 'high')) 
    & ((dgw.getisord_rf_month_february == 'low') | (dgw.getisord_rf_month_february == 'ns'))),'TC'] = 'TC'


dgs.loc[((dgs.getisord_trop_08 == 'high') & 
    ((dgs.getisord_cmaq_08 == 'low') | (dgs.getisord_cmaq_08 == 'ns')) 
    & ((dgs.getisord_rf_month_august == 'low') | (dgs.getisord_rf_month_august == 'ns'))),'TO'] = 'TO'

dgw.loc[((dgw.getisord_trop_02 == 'high') & 
    ((dgw.getisord_cmaq_02 == 'low') | (dgw.getisord_cmaq_02 == 'ns')) 
    & ((dgw.getisord_rf_month_february == 'low') | (dgw.getisord_rf_month_february == 'ns'))),'TO'] = 'TO'

dgs.loc[((dgs.getisord_rf_month_august == 'high') & 
    ((dgs.getisord_cmaq_08 == 'low') | (dgs.getisord_cmaq_08 == 'ns')) 
    & ((dgs.getisord_trop_08 == 'low') | (dgs.getisord_trop_08 == 'ns'))),'EO'] = 'EO'

dgw.loc[((dgw.getisord_rf_month_february == 'high') & 
    ((dgw.getisord_cmaq_02 == 'low') | (dgw.getisord_cmaq_02 == 'ns')) 
    & ((dgw.getisord_trop_02 == 'low') | (dgw.getisord_trop_02 == 'ns'))),'EO'] = 'EO'


#TO - Tropomi only
#EO - Eclipse only
#CL - CMAQ only
#TC - TropOMI + CMAQ
#SW - TropOMI + Eclipse
#CL2 - CMAQ + Eclipse
v = ['SW','TC','CL2','EO','TO','CL']
[print(len(dgw[dgw[v[i]] == v[i]]),len(dgs[dgs[v[i]] == v[i]])) for i in range(len(v))]

for i in range(len(v)):
    if len(dgw[dgw[v[i]] == v[i]]) <= 9:
        print('drop dgw '+v[i])
    if len(dgs[dgs[v[i]] == v[i]]) <= 9:
        print('drop dgs '+v[i])


# check if answers have < 5 gc
len(pd.merge(dgw[dgw.SW == 'SW'],dgs[dgs.SW == 'SW'],on = 'u'))
len(pd.merge(dgw[dgw.TC == 'TC'],dgs[dgs.TC == 'TC'],on = 'u'))
len(pd.merge(dgw[dgw.CL2 == 'CL2'],dgs[dgs.CL2 == 'CL2'],on = 'u'))


# check if answers have < 5 gc
len(pd.merge(dgw[dgw.EO == 'EO'],dgs[dgs.EO == 'EO'],on = 'u'))
len(pd.merge(dgw[dgw.TO == 'TO'],dgs[dgs.TO == 'TO'],on = 'u'))
len(pd.merge(dgw[dgw.CL == 'CL'],dgs[dgs.CL == 'CL'],on = 'u'))


#########################################################################################################
# Plot exploratory hotspots
# Plot clusters
# high confidence

fig,axs = plt.subplots(1,2,figsize=(8,5))
axs=axs.ravel()

#high
ax=axs[0]
dgs.plot(ax=ax,facecolor='lightgray')
dgs[dgs.HH == 'HH'].plot(ax=ax,facecolor='red',alpha=0.6)
#dgs[dgs.LL == 'LL'].plot(ax=ax,facecolor='mediumblue',alpha=0.6)
plt.rcParams['hatch.linewidth'] = 0.3
dgw[dgw.HH == 'HH'].plot(ax=ax,facecolor='red',alpha=0.6,hatch='//')
#dgw[dgw.LL == 'LL'].plot(ax=ax,facecolor='mediumblue',alpha=0.6,hatch='//')
class1.plot(ax=ax,edgecolor='k',alpha = 0.3)
ax.set_title('(a) High Agreement Hotspot')
ax.axis('off')

# low
axs=axs.ravel()
ax=axs[1]
dgs.plot(ax=ax,facecolor='lightgray')
dgs[dgs.SW == 'SW'].plot(ax=ax,facecolor='goldenrod',alpha=0.6,label=True)
dgs[dgs.CL == 'CL'].plot(ax=ax,facecolor='forestgreen',alpha=0.6)
#dgw[dgw.SW == 'SW'].plot(ax=ax,facecolor='green',alpha=0.6,hatch="//")
#cl2 = dgw[dgw.CL2 == 'CL2'][~dro]
#
f0 = pd.merge(dgs[dgs.CL == 'CL'],dgw[dgw.CL2 == 'CL2'], on = 'u')
f0['geometry'] = f0.geometry_x
#
dgw[dgw.CL2 == 'CL2'].plot(ax=ax,facecolor='green',alpha=0.3,hatch="//")
#f0[f0.CL2_y == 'CL2'].plot(ax=ax,facecolor='forestgreen',alpha=0.6,hatch="//")
#cl2.plot(ax=ax,facecolor='green',alpha=0.3,hatch="//")
plt.rcParams['hatch.linewidth'] = 0.3
class1.plot(ax=ax,edgecolor='k',alpha=0.3)
ax.set_title('(b) Robust Hotspots')
ax.axis('off')

plt.tight_layout()
plt.savefig('/projects/b1045/montgomery/microsoft/high_confidence_clusters_v2.png',dpi=300)
plt.show()



#
#########################################################################################################
######

fig,axs = plt.subplots(1,3,figsize=(12,5))
axs=axs.ravel()

#high
ax=axs[0]
dgs.plot(ax=ax,facecolor='lightgray')
dgs[dgs.HH == 'HH'].plot(ax=ax,facecolor='red',alpha=0.6)
#dgs[dgs.LL == 'LL'].plot(ax=ax,facecolor='mediumblue',alpha=0.6)
plt.rcParams['hatch.linewidth'] = 0.3
dgw[dgw.HH == 'HH'].plot(ax=ax,facecolor='red',alpha=0.6,hatch='//')
#dgw[dgw.LL == 'LL'].plot(ax=ax,facecolor='mediumblue',alpha=0.6,hatch='//')
class1.plot(ax=ax,edgecolor='k',alpha = 0.3)
ax.set_title('(a) Consensus Hotspot')
ax.axis('off')

# low
axs=axs.ravel()
ax=axs[1]
dgs.plot(ax=ax,facecolor='lightgray')
dgs[dgs.SW == 'SW'].plot(ax=ax,facecolor='goldenrod',alpha=0.6,label=True)
#dgs[dgs.CL2 == 'CL2'].plot(ax=ax,facecolor='green',alpha=0.6)
#dgs[dgs.TC == 'TC'].plot(ax=ax,facecolor='purple',alpha=0.6)

#dgw[dgw.SW == 'SW'].plot(ax=ax,facecolor='goldenrod',alpha=0.3,hatch="//")
dgw[dgw.CL2 == 'CL2'].plot(ax=ax,facecolor='green',alpha=0.3,hatch="//")
#dgw[dgw.TC == 'TC'].plot(ax=ax,facecolor='purple',alpha=0.3,hatch="//")


#cl2 = dgw[dgw.CL2 == 'CL2'][~dro]
cl2.plot(ax=ax,facecolor='green',alpha=0.6,hatch="//")
plt.rcParams['hatch.linewidth'] = 0.3
class1.plot(ax=ax,edgecolor='k',alpha=0.3)
ax.set_title('(b) Medium Agreement')
ax.axis('off')


# low
axs=axs.ravel()
ax=axs[2]
dgs.plot(ax=ax,facecolor='lightgray')
dgs[dgs.EO == 'EO'].plot(ax=ax,facecolor='darkorange',alpha=0.6,label=True)
dgs[dgs.CL == 'CL'].plot(ax=ax,facecolor='purple',alpha=0.6)
dgs[dgs.TO == 'TO'].plot(ax=ax,facecolor='blue',alpha=0.6)

#dgw[dgw.EO == 'EO'].plot(ax=ax,facecolor='darkorange',alpha=0.3,hatch="//")
#dgw[dgw.CL == 'CL'].plot(ax=ax,facecolor='purple',alpha=0.3,hatch="//")
dgw[dgw.TO == 'TO'].plot(ax=ax,facecolor='blue',alpha=0.3,hatch="//")

plt.rcParams['hatch.linewidth'] = 0.3
class1.plot(ax=ax,edgecolor='k',alpha=0.3)
ax.set_title('(c) Low Agreement')
ax.axis('off')


plt.tight_layout()
plt.savefig('/projects/b1045/montgomery/microsoft/high_confidence_clusters_v3_p01.png',dpi=300)
plt.show()

# how many people affected
dgs[~dgs['HH'].isna()]['total_pop'].sum()




# how many people affected
dgs[~dgs['HH'].isna()]['total_pop'].sum()



#
#########################################################################################################
# Plot exploratory hotspots
# Plot clusters
# high confidence
plt.rcParams['hatch.linewidth'] = 0.3

fig,axs = plt.subplots(2,2,figsize=(8,5))
axs=axs.ravel()
[dgs.plot(ax=axs[i],facecolor='lightgray') for i in range(4)]
[class1.plot(ax=axs[i],edgecolor='k',alpha = 0.3) for i in range(4)]
dgs[dgs.CL == 'CL'].plot(ax=axs[0],facecolor='forestgreen',alpha=0.6)
dgw[dgw.CL == 'CL'].plot(ax=axs[1],facecolor='green',alpha=0.6,hatch="//")
axs[0].set_title('August \n Modeled: CMAQ Only')
axs[1].set_title('February\n Modeled: CMAQ Only')

# low
dgs.plot(ax=ax,facecolor='lightgray')
dgs[dgs.CL2 == 'CL2'].plot(ax=axs[2],facecolor='forestgreen',alpha=0.6)
dgw[dgw.CL2 == 'CL2'].plot(ax=axs[3],facecolor='green',alpha=0.6,hatch="//")
axs[2].set_title('August \n Modeled: CMAQ + Eclipse')
axs[3].set_title('February\n Modeled: CMAQ + Eclipse')

[axs[i].axis('off') for i in range(4)]
plt.tight_layout()
#plt.savefig('/projects/b1045/montgomery/microsoft/updatedclusters.png',dpi=300)
plt.show()

# Check significance of other variables in hotspots
###########################################


cat_col = ['winddir_re', 'airt_rean', 'sf_rean', 'rh_rean', 'ws_rean','ndvi', 
        'arterials', 'highway', 'industry_a', 'res_ar', 'com_ar', 'manu_ar',
        'open_ar', 'oth_ar', 'avg_speed', 'speed_ct', 'avg_speed_',
        'total','winddir_rean02', 'airt_rean02', 'sf_rean02',
       'rh_rean02', 'ws_rean02',
       'total_NHL', 'White', 'Black', 'Native', 'Asian', 'HawaiiPI', 'Other',
       'Two', 'ALUKE010', 'ALUKE011', 'HispLatino','pWhite', 'pBlack', 'pHispLatino', 'pAsian',
        'incpc_all', 'pub_assist','ndvi_feb',
        'cmaq_08','trop_08','rf_month_august',
        'cmaq_02','trop_02','rf_month_february']


# Get average characteristics of each cluster
fout = pd.DataFrame([np.array(dgs[dgs.HH == 'HH'][cat_col].mean()),
np.array(dgs[(dgs.LL == 'LL')][cat_col].mean()),
np.array(dgs[(dgs.SW == 'SW')][cat_col].mean()),
np.array(dgs[(dgs.CL == 'CL')][cat_col].mean()),
np.array(dgs[cat_col].mean())]).T

wint=True
if wint:
    # Get average characteristics of each cluster
    fout = pd.DataFrame([np.array(dgw[dgw.HH == 'HH'][cat_col].mean()),
    np.array(dgw[(dgw.LL == 'LL')][cat_col].mean()),
    np.array(dgw[(dgw.SW == 'SW')][cat_col].mean()),
    np.array(dgw[(dgw.CL == 'CL')][cat_col].mean()),
    np.array(dgw[cat_col].mean())]).T


fout['chars'] = cat_col

fout.columns = ['High-High','Low-Low','SW Cluster','Loop Cluster','Chicago Avg.','Characteristics']
fout

df_HH = dgs[dgs['HH']=='HH']
df_LL = dgs[dgs['LL']=='LL']
df_CL = dgs[dgs['CL']=='CL']
df_SW = dgs[dgs['SW']=='SW']

# using t test to check sample
import scipy.stats as stats

#perform one sample t-test

stat_sig_HH = []; 
stat_sig_LL = []; 
stat_sig_CL = []; 
stat_sig_SW = []; 

def f_test(group1, group2):
    f = np.var(group1, ddof=1)/np.var(group2, ddof=1)
    nun = group1.size-1
    dun = group2.size-1
    p_value = 1-stats.f.cdf(f, nun, dun)
    return f, p_value


for x in cat_col:
    print(np.var(df_HH[x])**-1*np.var(dgs[x]))
    if f_test(df_HH[x],dgs[x])[1] < 0.05: boo = False
    else: boo = True
    stat_sig_HH.append(stats.ttest_ind(a=df_HH[x], b=dgs[x], equal_var=boo)[1])
    stat_sig_LL.append(stats.ttest_ind(a=df_LL[x], b=dgs[x], equal_var=boo)[1])
    stat_sig_CL.append(stats.ttest_ind(a=df_CL[x], b=dgs[x], equal_var=boo)[1])
    stat_sig_SW.append(stats.ttest_ind(a=df_SW[x], b=dgs[x], equal_var=boo)[1])

n = pd.DataFrame([cat_col,stat_sig_HH,stat_sig_LL,stat_sig_CL,stat_sig_SW]).T
n.columns = ['Characteristics','p_HH','p_LL','p_CL','p_SW']
nout = pd.merge(n,fout,on='Characteristics')

#nout['High-High'] = (nout['High-High'] - nout['Chicago Avg.'])/nout['Chicago Avg.']*100
#nout['SW Cluster'] = (nout['SW Cluster'] - nout['Chicago Avg.'])/nout['Chicago Avg.']*100
#nout['Loop Cluster'] = (nout['Loop Cluster'] - nout['Chicago Avg.'])/nout['Chicago Avg.']*100

cat = ['p_HH','p_LL','p_CL','p_SW']
cat2 = ['High-High','Low-Low','Loop Cluster','SW Cluster']

for i in range(len(cat)):
    
    tmp = nout[(nout[cat[i]] < 0.0125) &  (nout[cat[i]] >= 0.0)][cat2[i]].reset_index(drop=True)
    tmpfill = [str(round(tmp[i],2))+'*' for i in range(len(tmp))]
    nout.loc[(nout[cat[i]] < 0.0125) &  (nout[cat[i]] >= 0.0),cat2[i]] = tmpfill

#    tmp = nout[(nout[cat[i]] < 0.01) &  (nout[cat[i]] > 0.001)][cat2[i]].reset_index(drop=True)
#    tmpfill = [str(round(tmp[i],2))+'**' for i in range(len(tmp))]
#    nout.loc[(nout[cat[i]] < 0.01) &  (nout[cat[i]] > 0.001),cat2[i]] = tmpfill

#    tmp = nout[(nout[cat[i]] < 0.001) &  (nout[cat[i]] >= 0)][cat2[i]].reset_index(drop=True)
#    tmpfill = [str(round(tmp[i],2))+'***' for i in range(len(tmp))]
#    nout.loc[(nout[cat[i]] < 0.001) &  (nout[cat[i]] >= 0),cat2[i]] = tmpfill

pd.options.display.float_format = "{:.2f}".format

nout = nout.drop(['p_HH','p_LL','p_CL','p_SW','Low-Low'],axis=1)

droppy = ['manu_ar','open_ar','oth_ar','speed_ct','avg_speed_','pWhite','Other','ALUKE010','ALUKE011','total_NHL',
            'pBlack','pHispLatino','pNativeAmer','pAsian','HawaainPac','Other','Two+','sf_rean','rh_rean','winddir_re']

for i in range(len(droppy)):
    nout = nout[nout.Characteristics != droppy[i]].reset_index(drop=True)

nout.to_csv('/projects/b1045/montgomery/microsoft/sig_table_v2_full_wint.csv')



#########################################################################################################
# Supplementary Figures
#########################################################################################################
# Plot raw data
# 
plt.rc('legend', fontsize=8) 
import matplotlib
matplotlib.rcParams.update({'font.size': 8})

fig,axs = plt.subplots(2,3,figsize=(7,5))
axs=axs.ravel()

x,y = dg.geometry[0].centroid.x-.27,dg.geometry[0].centroid.y-0.02

vmins = [5,5,5,5,10,5]
vmaxs = [20,10,15,20,25,20]
units=['ppb',r'VCD','ppb']*2

t = ['cmaq_08','trop_08','rf_month_august','cmaq_02','trop_02','rf_month_february']
for i in range(len(t)):
    txt = r"$\mu$ = %.1f"%(dg[t[i]].median())+"\n"+r"$\sigma$ = %.1f"%(dg[t[i]].std())
    axs[i].text(s=txt, x = x,y=y,fontsize=9,zorder=100)
    dg.plot(t[i],ax=axs[i],legend=True,vmin = vmins[i],vmax=vmaxs[i],
        legend_kwds={"shrink":.5, 'ticks': [vmins[i],(vmins[i]+vmaxs[i])/2,vmaxs[i]],'label':units[i]})
    axs[i].axis('off')

fig.text(.33/2, 0.98, 'CMAQ', ha='center',fontsize=13)
fig.text((0.66+0.33)/2, 0.98, 'TropOMI', ha='center',fontsize=13)
fig.text((1+0.66)/2, 0.98, 'Eclipse', ha='center',fontsize=13)

fig.text(0.01, 0.25, 'February 2022', va='center', rotation='vertical',fontsize=13)
fig.text(0.01, 0.75, 'August 2021', va='center', rotation='vertical',fontsize=13)

plt.tight_layout()
plt.savefig('Eclipse-CMAQ-Trop-Chi.png',dpi=300)
plt.show()

#########################################################################################################
# Plot raw data 1:1 and lin reg

from sklearn.linear_model import LinearRegression
def get_linreg(x,y):
    model = LinearRegression()
    x,y = np.array(x),np.array(y)
    x = x.reshape(-1,1) #np.array([5, 15, 25, 35, 45, 55]).reshape((-1, 1))
    y = y.reshape(-1) #np.array([5, 20, 14, 32, 22, 38])3
    model = LinearRegression().fit(x, y)
    r_sq = model.score(x, y)
    print(f"intercept: {model.intercept_}")
    y_pred = model.predict(np.arange(0,35).reshape(-1,1))
    return y_pred

fig,axs = plt.subplots(1,3,figsize=(8,3))

ti = ['August 2021', 'February 2022']
t = [['cmaq_08','trop_08','rf_month_august'],['cmaq_02','trop_02','rf_month_february']]
markerstyle = ['o']+['^']
c = ['navy','darkgreen']
mo = ['Aug.','Feb.']

for i in range(len(t)):
    axs[0].scatter(dg[t[i][1]],dg[t[i][0]],alpha = 0.5,s=10,c=c[i],marker=markerstyle[i])
    axs[0].plot([0,35],[0,35],c='k',zorder=0.1,alpha=0.7)
    axs[0].set_xlim([0,35]);axs[0].set_ylim([0,35]);
    axs[0].set_box_aspect(1)
    #axs[0].set_title('1:1 for %s'%(ti[i]))
    axs[0].set_ylabel('CMAQ')
    axs[0].set_xlabel('TropOMI')
    y_pred = get_linreg(np.array(dg[t[i][1]]),np.array(dg[t[i][0]]))
    a = np.corrcoef(dg[t[i][1]],dg[t[i][0]])[0][1]
    b = np.mean(dg[t[i][1]]-dg[t[i][0]])*-1
    axs[0].plot(np.arange(0,35),y_pred,zorder=0.1,c=c[i],linestyle='--',alpha=0.5,label='%s: r = %.1f, mb = %.1f'%(mo[i],a,b))
    #
    #axs[0].annotate(xy=(3,1),text="r = %.1f, mb = %.1f"%(a,b))
    #
    axs[1].scatter(dg[t[i][2]],dg[t[i][1]],alpha = 0.5,s=10,c=c[i],marker=markerstyle[i])
    axs[1].plot([0,35],[0,35],c='k',zorder=.1,alpha=0.7)
    axs[1].set_xlim([0,35]);axs[1].set_ylim([0,35]);
    axs[1].set_box_aspect(1)
    #axs[1].set_title('1:1 for %s'%(ti[i]))
    axs[1].set_ylabel('TropOMI')
    axs[1].set_xlabel('Eclipse')
    #
    a = np.corrcoef(dg[t[i][2]],dg[t[i][1]])[0][1]
    b = np.mean(dg[t[i][2]]-dg[t[i][1]])*-1
    y_pred = get_linreg(dg[t[i][2]],dg[t[i][1]])
    axs[1].plot(np.arange(0,35),y_pred,zorder=0.1,c=c[i],linestyle='--',alpha=0.5,label='%s: r = %.1f, mb = %.1f'%(mo[i],a,b))
    #axs[1].annotate(xy=(3,1),text="r = %.1f, mb = %.1f"%(a,b))
    #
    #
    axs[2].scatter(dg[t[i][2]],dg[t[i][0]],alpha = 0.5,s=10,c=c[i],marker=markerstyle[i])
    axs[2].plot([0,35],[0,35],c='k',zorder=.1,alpha=0.7)
    axs[2].set_xlim([0,35]);axs[2].set_ylim([0,35]);
    axs[2].set_box_aspect(1)
    #axs[2].set_title('1:1 for %s'%(ti[i]))
    axs[2].set_ylabel('CMAQ')
    axs[2].set_xlabel('Eclipse')
    #
    a = np.corrcoef(dg[t[i][2]],dg[t[i][0]])[0][1]
    b = np.mean(dg[t[i][2]]-dg[t[i][0]])*-1
    y_pred = get_linreg(dg[t[i][2]],dg[t[i][0]])
    axs[2].plot(np.arange(0,35),y_pred,zorder=0.1,c=c[i],linestyle='--',alpha=0.5,label='%s: r = %.1f, mb = %.1f'%(mo[i],a,b))
    #
    #
    #axs[2].annotate(xy=(3,1),text="r = %.1f, mb = %.1f"%(a,b))
    

[axs[i].legend(loc='upper center', bbox_to_anchor=(0.5, 1.23),
          ncol=1, fancybox=True, prop={'size': 8}) for i in range(3)]
plt.tight_layout()
plt.savefig('Eclipse-CMAQ-Trop-1.png',dpi=300)
plt.show()
#high

#########################################################################################################
# Do daily clustering
# add daily eclipse
d = '/projects/b1045/montgomery/microsoft/'
daily = pd.read_csv('/projects/b1045/montgomery/microsoft/daily_data_grid_eclipse_rf.csv')
tmp = gpd.GeoDataFrame([dg.u,dg.geometry],crs = dg.crs).T
daily = gpd.GeoDataFrame(pd.merge(tmp,daily,on='u'),crs = dg.crs)
daily['dt'] = pd.to_datetime(daily.date).dt.date


t = daily.date.unique()
ts = t[:31]
tw = t[31:]
trop_08_daily = pd.read_csv(d+'trop_daily_08.csv',index_col=0)
trop_02_daily = pd.read_csv(d+'trop_daily_02.csv',index_col=0)
trop_08_daily['date'] = np.array([[ts[i]]*986 for i in range(31)]).ravel()
trop_02_daily['date'] = np.array([[tw[i]]*986 for i in range(28)]).ravel()
trop = pd.concat([trop_08_daily,trop_02_daily]).reset_index(drop=True)
trop.columns = ['u','d','trop','date']

din = pd.read_csv(d+'eclipse_daily_1pm.csv',index_col=0)
din['dt'] = pd.to_datetime(din.date).dt.date
daily = gpd.GeoDataFrame(pd.merge(din,daily,on=['u','dt']),crs = dg.crs)


#add daily cmaq
from netCDF4 import Dataset
no2_08 = np.array([np.mean(Dataset(d+'all_no2.nc')['NO2'][24*i:(24*i+24)],axis=0)[0] for i in range(30)]+[np.mean(Dataset(d+'all_no2.nc')['NO2'][24*30:],axis=0)[0]])
no2_02 = np.array([np.mean(Dataset(d+'all_no2_202202.nc')['NO2'][24*i:(24*i+24)],axis=0)[0] for i in range(28)])

# 1 PM ct
x,y = 18,2
no2_08 = np.array([np.mean(Dataset(d+'all_no2.nc')['NO2'][(24*i+x):(24*i+x+y)],axis=0)[0] for i in range(30)]+[np.mean(Dataset(d+'all_no2.nc')['NO2'][24*30:],axis=0)[0]])
no2_02 = np.array([np.mean(Dataset(d+'all_no2_202202.nc')['NO2'][(24*i+x+y):(24*i+y)],axis=0)[0] for i in range(28)])

# 1 PM ct -- COLUMN
x,y = 18,2
di = '/projects/b1045/wrf-cmaq/output/Chicago_LADCO/output_202108_1.33km_sf_rrtmg_5_8_1_v3852.bk/column/'
no2_08_col = np.array([np.mean(Dataset(di+'no2_column_202108.nc')['NO2'][(24*i+x):(24*i+x+y)],axis=0)[0] for i in range(30)]+[np.mean(Dataset(d+'all_no2.nc')['NO2'][24*30:],axis=0)[0]])
no2_02_col = np.array([np.mean(Dataset(d+'all_no2_202202.nc')['NO2'][(24*i+x+y):(24*i+y)],axis=0)[0] for i in range(28)])


t = daily.date_y.unique()

lonuf=[98, 132] #indices outside of Chicago
latuf=[123, 152] # indices of other corner outside of chicago

no2_08_crop = [np.array(no2_08[i][lonuf[0]:lonuf[1]]) for i in range(31)]
no2_08_crop = np.array([np.array(no2_08_crop[i].T[latuf[0]:latuf[1]].T) for i in range(31)])
no2_08_crop = np.array([no2_08_crop[t].ravel() for t in range(31)])

no2_08col_crop = [np.array(no2_08_col[i][lonuf[0]:lonuf[1]]) for i in range(31)]
no2_08col_crop = np.array([np.array(no2_08col_crop[i].T[latuf[0]:latuf[1]].T) for i in range(31)])
no2_08col_crop = np.array([no2_08col_crop[t].ravel() for t in range(31)])


no2_02_crop = [np.array(no2_02[i][lonuf[0]:lonuf[1]]) for i in range(28)]
no2_02_crop = np.array([np.array(no2_02_crop[i].T[latuf[0]:latuf[1]].T) for i in range(28)])
no2_02_crop = np.array([no2_02_crop[t].ravel() for t in range(28)])

tmp_shp = gpd.read_file(d+'rectangle_chicago_1km_grid.shp')

f = pd.DataFrame()
for i in range(31):
    t2 = gpd.GeoDataFrame(pd.DataFrame([tmp_shp.u,tmp_shp.geometry]).T,crs=tmp_shp.crs)
    t2['no2_cmaq_crop'] = no2_08_crop[i]
    t2['no2_cmaq_col'] = no2_08col_crop[i]
    t2['date'] = ts[i]
    f=f.append(t2)


for i in range(28):
    t2 = gpd.GeoDataFrame(pd.DataFrame([tmp_shp.u,tmp_shp.geometry]).T,crs=tmp_shp.crs)
    t2['no2_cmaq_crop'] = no2_02_crop[i]
    t2['date'] = t[i+31]
    f=f.append(t2)

daily['date'] = daily.date_y
daily = pd.merge(daily,f,on=['u','date'])
daily['geometry'] = daily.geometry_x
daily = daily.drop(['geometry_x','geometry_y'],axis=1)
daily = pd.merge(daily,trop, on = ['u','date'])
daily = gpd.GeoDataFrame(daily,crs=dg.crs)
#daily['trop'] = trop['2']

#start daily clustering
w_surface = weights.Queen.from_dataframe(dg)

t0 = t[:31]
final = pd.DataFrame()
for i in range(len(t0)):
    tmp = daily[daily.date_y == t0[i]].reset_index(drop=True)
    final = final.append(get_ordg_hotspots_from_shp2(tmp,w_surface,['cal_no2','no2_cmaq_crop','no2_cmaq_col','trop']))

finals = gpd.GeoDataFrame(final,crs=dg.crs)

t0 = t[31:]
final = pd.DataFrame()
for i in range(len(t0)):
    tmp = daily[daily.date == t0[i]].reset_index(drop=True)
    final = final.append(get_ordg_hotspots_from_shp2(tmp,w_surface,['rf_no2','no2_cmaq_crop','trop']))

finalw = gpd.GeoDataFrame(final,crs=dg.crs)

monthly = daily.groupby(['u']).mean().reset_index()
monthly = gpd.GeoDataFrame(pd.merge(monthly,tmp,on='u'),crs = tmp.crs)
monthly = get_ordg_hotspots_from_shp2(monthly,w_surface,['cal_no2','no2_cmaq_crop','no2_cmaq_col','trop'])



#==  analyze successful retrievals over moth
p = ['getisord_cal_no2','getisord_no2_cmaq_crop','getisord_trop','CO']
z = ['(a) Eclipse','(b) WRF-CMAQ Surface','(c) TropOMI','(d) Consensus']

finals['CO'] = ''
finals.loc[(finals[p[0]]=='high') &  (finals[p[1]] == 'high') &  (finals[p[2]] == 'high'),'CO'] = 'high'

monthly['CO'] = ''
monthly.loc[(monthly[p[0]]=='high') &  (monthly[p[1]] == 'high') &  (monthly[p[2]] == 'high'),'CO'] = 'high'

fig,ax = plt.subplots(1,4,figsize=(6,3))
i = 0
for j in range(len(p)):
                crop = monthly
                ch = crop[crop[p[j]] == 'high']
        #        cns = crop[crop[p[j]] == 'ns']
                crop.plot(facecolor='gray',ax=ax[j])
                ch.plot(facecolor='red',ax=ax[j])
                ax[j].axis('off')
                ax[j].set_title(str(z[j]), fontsize=8)
                

plt.savefig('daily_1pm_082021_surface_all3_onemonth.png',dpi=350)
plt.tight_layout()
plt.show()


#==  analyze successful retrievals over column for month
p = ['getisord_no2_cmaq_col','getisord_trop','CO']
z = ['(a) WRF-CMAQ Surface','(b) TropOMI','(c) Consensus']

monthly['CO'] = ''
monthly.loc[(monthly[p[0]]=='high') &  (monthly[p[1]] == 'high') &  (monthly[p[1]] == 'high'),'CO'] = 'high'

fig,ax = plt.subplots(1,3,figsize=(6,3))
i = 0
for j in range(len(p)):
                crop = monthly
                ch = crop[crop[p[j]] == 'high']
        #        cns = crop[crop[p[j]] == 'ns']
                crop.plot(facecolor='gray',ax=ax[j])
                ch.plot(facecolor='red',ax=ax[j])
                ax[j].axis('off')
                ax[j].set_title(str(z[j]), fontsize=8)
                

plt.savefig('daily_1pm_082021_column_all3_onemonth.png',dpi=350)
plt.tight_layout()
plt.show()


p = ['getisord_cal_no2','getisord_no2_cmaq_crop','getisord_trop','CO']
z = ['Eclipse','WRF-CMAQ Surface','TropOMI','Consensus']

finals['CO'] = ''
finals.loc[(finals[p[0]]=='high') &  (finals[p[1]] == 'high') &  (finals[p[2]] == 'high'),'CO'] = 'high'

fig,axs = plt.subplots(12,4,figsize=(6,10))
i = 0
for t in range(len(ts)):
    crop = finals[finals.date == ts[t]].reset_index()
    if len(crop) > 0:
            ax = axs[i]
            for j in range(len(p)):
                ch = crop[crop[p[j]] == 'high']
        #        cns = crop[crop[p[j]] == 'ns']
                crop.plot(facecolor='gray',ax=ax[j])
                ch.plot(facecolor='red',ax=ax[j])
                ax[j].axis('off')
                if t == 0: ax[j].set_title(str(z[j]) + '\n'+ str(crop.dt[0]), fontsize=6)
                else: ax[j].set_title(str(crop.dt[0]), fontsize=6)
#
            i = i+1

plt.savefig('daily_1pm_082021_surface_all3.png',dpi=350)
plt.tight_layout()
plt.show()


p = ['getisord_no2_cmaq_col','getisord_trop','CO']
z = ['WRF-CMAQ Column','TropOMI','WRF-CMAQ Column U TropOMI']

finals['CO'] = ''
finals.loc[(finals[p[0]]=='high') &  (finals[p[1]] == 'high') &  (finals[p[1]] == 'high'),'CO'] = 'high'

fig,axs = plt.subplots(12,3,figsize=(6,10))
i = 0
for t in range(len(ts)):
    crop = finals[finals.date == ts[t]].reset_index()
    if len(crop) > 0:
            ax = axs[i]
            for j in range(len(p)):
                ch = crop[crop[p[j]] == 'high']
        #        cns = crop[crop[p[j]] == 'ns']
                crop.plot(facecolor='gray',ax=ax[j])
                ch.plot(facecolor='red',ax=ax[j])
                ax[j].axis('off')
                if t == 0: ax[j].set_title(str(z[j]) + '\n'+ str(crop.dt[0]), fontsize=6)
                else: ax[j].set_title(str(crop.dt[0]), fontsize=8)
#
            i = i+1

plt.savefig('daily_1pm_082021_column.png',dpi=350)
plt.tight_layout()
plt.show()


def get_corr(a,b):
    a = np.array(a)
    b = np.array(b)
    b = b[~np.isnan(a)]
    a = a[~np.isnan(a)]
    return np.corrcoef(a,b)[0][1]

get_corr(finals.trop,finals.no2_cmaq_col)
get_corr(finals.trop,finals.no2_cmaq_crop)
mo = finals.groupby(['u']).mean().reset_index()
get_corr(mo.trop,mo.no2_cmaq_col)
get_corr(mo.trop,mo.no2_cmaq_crop)


# get intersection
d = '/projects/b1045/montgomery/microsoft/'
fin = pd.read_csv(d+'ec_on_cmaq_wint.csv',index_col=0)
fin = fin.groupby('hour').mean()
daily_winddir = np.array([fin.winddir_rean[24*i:(24*i+24)].mean() for i in range(30)] + [fin.winddir_rean[24*30:(24*30+24)].mean()])

t8 = t[0:31]
t2 = t[31:]

d3 = t8[(daily_winddir < 360) & (daily_winddir >=270)]
d2 = t8[(daily_winddir < 270) & (daily_winddir >=180)]
d1 = t8[(daily_winddir < 180) & (daily_winddir >=90)]
d0 = t8[(daily_winddir < 90) & (daily_winddir >=0)]

import pandas as pd
import numpy as np
from netCDF4 import Dataset
#import wrf
# get intersection
d = '/projects/b1045/montgomery/microsoft/th_2022.nc'
fin = Dataset(d)
xx,yy = np.meshgrid(fin['lon'],fin['lat'])
m = mask_given_shapefile(xx,yy,dgs)
daily_wind = np.array([fin['wind_from_direction'][t][m].mean() for t in range(32,60)])

pd.DataFrame(daily_wind).to_csv('/projects/b1045/montgomery/microsoft/daily_winddir_feb_chicago.csv')

d3w = t2[(daily_wind < 360) & (daily_wind >=270)]
d2w = t2[(daily_wind < 270) & (daily_wind >=180)]
d1w = t2[(daily_wind < 180) & (daily_wind >=90)]
d0w = t2[(daily_wind < 90) & (daily_wind >=0)]

def make_wind_df(d0,final,tmp_shp,v,w_surface):
    final_d0 = pd.DataFrame()
    for i in range(len(d0)):
        final_d0 = final_d0.append(final[final.date == d0[i]])
    #
    final_d0 = final_d0.groupby('u').mean().reset_index()
    final_d0 = gpd.GeoDataFrame(pd.merge(tmp_shp,final_d0,on='u'),crs=tmp_shp.crs)
    # add hotspots
    final_d0 = get_ordg_hotspots_from_shp(final_d0,w_surface,v)
    #find intersection
    final_d0['HH'] = ''
    final_d0['LL'] = ''
    final_d0.loc[(final_d0.getisord_rf_no2 == 'high') & (final_d0.getisord_no2_cmaq_crop == 'high'),'HH'] = 'HH'
    final_d0.loc[(final_d0.getisord_rf_no2 == 'low') & (final_d0.getisord_no2_cmaq_crop == 'low'),'LL'] = 'LL'
    #
    return final_d0

v= ['rf_no2','no2_cmaq_crop']
final_d0 = make_wind_df(d0,finals,tmp_shp,v,w_surface)
final_d1 = make_wind_df(d1,finals,tmp_shp,v,w_surface)
final_d2 = make_wind_df(d2,finals,tmp_shp,v,w_surface)
final_d3 = make_wind_df(d3,finals,tmp_shp,v,w_surface)

final_d0w = make_wind_df(d0w,finalw,tmp_shp,v,w_surface)
final_d1w = make_wind_df(d1w,finalw,tmp_shp,v,w_surface)
final_d2w = make_wind_df(d2w,finalw,tmp_shp,v,w_surface)
final_d3w = make_wind_df(d3w,finalw,tmp_shp,v,w_surface)


# plot the wind direfction clusters
x,y = dg.geometry[0].centroid.x-.28,dg.geometry[0].centroid.y-0.05
##########################################################################################

HH08_extent = gpd.GeoDataFrame({'geometry':gpd.GeoSeries(unary_union(dgs[dgs.HH == 'HH'].geometry))},crs=dgs.crs)
HH02_extent = gpd.GeoDataFrame({'geometry':gpd.GeoSeries(unary_union(dgw[dgw.HH == 'HH'].geometry))},crs=dgs.crs)

#CL208_extent = gpd.GeoDataFrame({'geometry':gpd.GeoSeries(unary_union(dgs[dgs.CL2 == 'CL2'].geometry))},crs=dgs.crs)
CL202_extent = gpd.GeoDataFrame({'geometry':gpd.GeoSeries(unary_union(dgw[dgw.CL2 == 'CL2'].geometry))},crs=dgs.crs)

datas = [final_d0,final_d1,final_d2,final_d3]
datas = [gpd.clip(datas[i],chi_outer) for i in range(len(datas))]

fig,ax=plt.subplots(2,2,figsize=(4,4))
ax=ax.ravel()
texts = ['(a)\nn = %s\n0 <= WD < 90'%(str(len(d0))),
        '(b)\nn = %s\n90 <= WD < 180'%(str(len(d1))),
        '(c)\nn = %s\n180 <= WD < 270'%(str(len(d2))),
        '(d)\nn = %s\n270 <= WD < 360'%(str(len(d3)))]


for i in range(len(datas)):
    datas[i].plot(facecolor='lightgray',ax=ax[i])
    hi = datas[i][datas[i].HH == 'HH']
    #lo = datas[i][datas[i].LL == 'LL']
    if len(hi) > 0: hi.plot(facecolor='red',ax=ax[i])
    #if len(lo) > 0: lo.plot(facecolor='blue',ax=ax[i])
    ax[i].axis('off')
    #class1.plot(ax=ax[i],edgecolor='k',linewidth=0.9,alpha=0.3)
    HH08_extent.plot(ax=ax[i],edgecolor='k',linewidth=0.9,alpha=0.5,facecolor="None")
    ax[i].text(s=texts[i], x = x,y=y,fontsize=9,zorder=100)

plt.savefig('cluster_by_winddir_08.png',dpi=300)
plt.show()

##########################################################################################
datas = [final_d0w,final_d1w,final_d2w,final_d3w]
datas = [gpd.clip(datas[i],chi_outer) for i in range(len(datas))]

fig,ax=plt.subplots(2,2,figsize=(4,4))
ax=ax.ravel()
texts = ['(a)\nn = %s\n0 <= WD < 90'%(str(len(d0w))),
        '(b)\nn = %s\n90 <= WD < 180'%(str(len(d1w))),
        '(c)\nn = %s\n180 <= WD < 270'%(str(len(d2w))),
        '(d)\nn = %s\n270 <= WD < 360'%(str(len(d3w)))]

for i in range(len(datas)):
    datas[i].plot(facecolor='lightgray',ax=ax[i])
    hi = datas[i][datas[i].HH == 'HH']
    #lo = datas[i][datas[i].LL == 'LL']
    if len(hi) > 0: hi.plot(facecolor='red',ax=ax[i])
    #if len(lo) > 0: lo.plot(facecolor='blue',ax=ax[i])
    ax[i].axis('off')
    class1.plot(ax=ax[i],edgecolor='k',linewidth=0.9,alpha=0.3)
    HH02_extent.plot(ax=ax[i],edgecolor='k',linewidth=0.9,alpha=0.5,facecolor="None")
    CL202_extent.plot(ax=ax[i],edgecolor='forestgreen',linewidth=0.9,alpha=0.6,facecolor="None")
    ax[i].text(s=texts[i], x = x,y=y,fontsize=9,zorder=100)

plt.savefig('cluster_by_winddir_02.png',dpi=300)
plt.show()

##########################################################################################

# Check tropOMI counts
##########################
dg['trop_count_08v2'] = dg.trop_count_08/31*100
dg['trop_count_02v2'] = dg.trop_count_02/28*100

# plot valid trop pixels
fig,ax = plt.subplots(1,2,figsize=(7,4))
dg.plot('trop_count_08v2',vmin=0,vmax=100,ax=ax[0],legend=True,legend_kwds={'shrink':0.5})
ax[0].set_title('TropOMI: August 2021 \n Average Pixel Coverage = %.1f'%(np.mean(dg.trop_count_08v2))+'%')
dg.plot('trop_count_02v2',vmin=0,vmax=100,ax=ax[1],legend=True,legend_kwds={'shrink':0.5})
ax[1].set_title('TropOMI: February 2022 \n Average Pixel Coverage = %.1f'%(np.mean(dg.trop_count_02v2))+'%')
ax[0].axis('off');ax[1].axis('off');
plt.tight_layout()
plt.savefig('tropomi_valid_count.png',dpi=300)
plt.show()




# Plot station locations over Chicago
########################################

from cartopy.io.img_tiles import Stamen
from cartopy import crs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

streets = gpd.read_file('/projects/b1045/montgomery/shapefiles/geo_export_28abe032-5cb4-4d93-9c9e-cf51d26e8169.shp')
class1 = streets[streets['class'] == '1'].reset_index(drop=True)
chi_shp = gpd.read_file('/projects/b1045/montgomery/shapefiles/Chicago/commareas/geo_export_77af1a6a-f8ec-47f4-977c-40956cd94f97.shp')
eclipse = pd.read_csv('ec_on_cmaq.csv')
ec = eclipse[~eclipse['ec_NO2'].isna()].groupby(['lon','lat']).mean().reset_index()

fig, ax = plt.subplots(subplot_kw={'projection': crs.PlateCarree()})
tiler = Stamen(style='toner-lite')
ax.add_image(tiler, 11)
#ax.coastlines('10m')
#ax.add_feature(cartopy.cfeature.LAND,facecolor='None',edgecolor='k')
chi_shp.plot(edgecolor='k',facecolor='None',linewidth=0.5,ax=ax,alpha=0.7)
ax.scatter(ec.lon,ec.lat,c = 'k',s=10)
class1.plot(ax=ax,edgecolor='navy',linewidth=0.9)
ax.set_title('Eclipse Stations in Chicago')
#plt.show()

# add label
ax.text(0.05,0.82,s="O'Hare\n   ✈",transform=ax.transAxes,fontsize = 10,alpha=0.8)
ax.text(0.33,0.4,s='Midway\n   ✈',transform=ax.transAxes,fontsize = 10,alpha=0.8)
ax.text(0.4,0.05,s='EPA ComED',transform=ax.transAxes,fontsize = 10,alpha=0.8)

ax.text(0.82,0.305,s="I-90",transform=ax.transAxes,fontsize = 9,alpha=0.8, style='italic',c='navy')
ax.text(0.535,0.605,s='I-290',transform=ax.transAxes,fontsize = 9,alpha=0.8, style='italic',c='navy')
ax.text(0.54,0.45,s='I-55',transform=ax.transAxes,fontsize = 9,alpha=0.8, style='italic',c='navy')
ax.text(0.49,0.81,s='I-94',transform=ax.transAxes,fontsize = 9,alpha=0.8, style='italic',c='navy')
ax.text(0.66,0.33,s='I-94',transform=ax.transAxes,fontsize = 9,alpha=0.8, style='italic',c='navy')
ax.text(0.59,0.14,s='I-57',transform=ax.transAxes,fontsize = 9,alpha=0.8, style='italic',c='navy')
ax.text(0.8,0.14,s='I-94',transform=ax.transAxes,fontsize = 9,alpha=0.8, style='italic',c='navy')
ax.text(0.8,0.5,s='LSD',transform=ax.transAxes,fontsize = 9,alpha=0.8, style='italic',c='navy')
ax.text(0.73,0.71,s='LSD',transform=ax.transAxes,fontsize = 9,alpha=0.8, style='italic',c='navy')

ax.text(0.8,0.8,s='  Lake\nMichigan',transform=ax.transAxes,fontsize = 12,alpha=0.5, style='italic')

gl = ax.gridlines(crs=crs.PlateCarree(), draw_labels=True)
gl.xlabels_top = False
gl.ylabels_left = False
gl.xlines = False
gl.ylines = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

gl.ylocator = mticker.FixedLocator([41.6, 41.7, 41.8, 41.9, 42.0])
gl.xlocator = mticker.FixedLocator([-87.9, -87.8,-87.7, -87.6, -87.5])

#plt.savefig('/projects/b1045/montgomery/paper1/counties_with_95_population.pdf')
plt.savefig('NO2_prelim_stnlocation.png',dpi=300)
plt.show()


# Plot station locations over Chicago
######################################################################

def find_index_1d(stn_lon, stn_lat, wrf_lon, wrf_lat):
#stn_lon, stn_lat, wrf_lon, wrf_lat=    fin.lon, fin.lat, eclipse.lon, eclipse.lat
# stn -- points 
# wrf -- list
#for iz in range(1):
    xx=[];yy=[]
    for i in range(len(stn_lat)):
    #for i in range(1):
        abslat = np.abs(wrf_lat-stn_lat[i])
        abslon= np.abs(wrf_lon-stn_lon[i])
        c = np.maximum(abslon,abslat)
        latlon_idx = np.argmin(c)
        x = np.where(c == np.min(c))
        #add indices of nearest wrf point station
        xx.append(x[0][0])
        #yy.append(y)
    #
    #xx=[xx[i][0] for i in range(len(xx))];yy=[yy[i][0] for i in range(len(yy))]
    #return indices list
    return xx#, yy

distx = find_index_1d(eclipse.lon, eclipse.lat,fin.lon, fin.lat)
fin['u'][distx]
gin = gpd.read_file(d+'all_data_grid_eclipse_idw_idp1.geojson')
gin = pd.read_file(d+'')
fed = gin[gin.u.isin(fin['u'][distx])]


distx,disty = np.array(distx),np.array(disty)
c = np.max(distx,disty)

f = ckdnearest(eclipse,ok,gdfB_cols=['AADT'])

def ckdnearest(gdfA, gdfB, gdfB_cols=['Place']):
    A = np.concatenate(
        [np.array(geom.coords) for geom in gdfA.geometry.to_list()])
    B = [np.array(geom.coords) for geom in gdfB.geometry.to_list()]
    B_ix = tuple(itertools.chain.from_iterable(
        [itertools.repeat(i, x) for i, x in enumerate(list(map(len, B)))]))
    B = np.concatenate(B)
    ckd_tree = cKDTree(B)
    dist, idx = ckd_tree.query(A, k=1)
    idx = itemgetter(*idx)(B_ix)
    gdf = pd.concat(
        [gdfA, gdfB.loc[idx, gdfB_cols].reset_index(drop=True),
         pd.Series(dist, name='dist')], axis=1)
    return gdf

