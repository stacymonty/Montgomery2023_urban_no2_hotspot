microsoft_check.py

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
from scipy.stats import pearsonr
import cartopy
from cartopy import crs
from cartopy.io.img_tiles import Stamen
from netCDF4 import Dataset
from matplotlib.dates import DateFormatter
import matplotlib.dates as mdates
import geopandas as gpd
import cartopy.feature as cfeature


# functions
def stats_normalized(data,prediction):
    from sklearn.metrics import mean_squared_error as mse
    x,y=data[~np.isnan(data)],prediction[~np.isnan(data)] # get rid of NaNs
    x,y=data[~np.isnan(prediction)],prediction[~np.isnan(prediction)] # get rid of NaNs
    mu_d,mu_p = np.mean(x),np.mean(y)
    mb = np.mean(y-x)
    nmb = np.sum(y-x)/np.sum(x)*100
    rmse = np.sqrt(mse(x,y))
    nrmse = np.sqrt(mse(x,y)/np.mean(x))
    r,p = st.pearsonr(x,y)
    return mu_d,mu_p,mb,nmb,rmse,nrmse,r


f = pd.read_csv('/projects/b1045/montgomery/CMAQcheck/NO2_d03_2021_8_EPA_CMAQ_Combine.csv',index_col=0)
f = f[f.x != 0].reset_index(drop=True)
f['dt'] = pd.to_datetime(f['dt'])
f['h'] = [f.dt[i].hour for i in range(len(f))]
#f = f[(f['h']>=2) & (f['h']<=7)].reset_index(drop=True) # get different time slices
f = f[~pd.isna(f['Sample Measurement'])].reset_index(drop=True)
x =  np.array(f['Sample Measurement'])
y =  np.array(f['CMAQ'])
sdc = len(f['Site Num'].unique())
summ = stats_normalized(x,y) # returns stats for that 1 file


fc = f[(f['Site Num'] == 76) | (f['Site Num'] == 3103) | (f['Site Num'] == 219)].reset_index(drop=True)
fc = fc[~pd.isna(f['Sample Measurement'])].reset_index(drop=True)
x =  np.array(fc['Sample Measurement'])
y =  np.array(fc['CMAQ'])
snc = len(fc['Site Num'].unique())
csum = stats_normalized(x,y)

f = pd.read_csv('/projects/b1045/montgomery/CMAQcheck/NO2_d03_2022_2_EPA_CMAQ_Combine.csv',index_col=0)
f = f[f.x != 0].reset_index(drop=True)
f = f[~pd.isna(f['Sample Measurement'])].reset_index(drop=True)
f['dt'] = pd.to_datetime(f['dt'])
f['h'] = [f.dt[i].hour for i in range(len(f))]
f = f[(f['h']>=2) & (f['h']<=7)].reset_index()
x =  np.array(f['Sample Measurement'])
y =  np.array(f['CMAQ'])
wds = len(f['Site Num'].unique())
wint = stats_normalized(x,y)

fc = f[(f['Site Num'] == 76) | (f['Site Num'] == 3103) | (f['Site Num'] == 219)].reset_index(drop=True)
fc = fc[~pd.isna(f['Sample Measurement'])].reset_index(drop=True)
x =  np.array(fc['Sample Measurement'])
y =  np.array(fc['CMAQ'])
wnc = len(fc['Site Num'].unique())
cwint = stats_normalized(x,y)


stats = pd.DataFrame({'sum':summ,'chi_sum':csum,'wint':wint,'chi_wint':cwint})
stats = stats.T
stats.columns = ['mu_d','mu_p','mb','nmb','rmse','nrmse','r']
stats['n'] = [sdc,snc,wds,wnc]
stats.to_csv('microsoft_check.csv')


lau,lonu = f.Latitude.unique(),f.Longitude.unique()

fig,ax = plt.subplots(4,2,figsize=(9,9),sharex=True)
ax = ax.ravel()
plt.xticks(rotation=90)
stats = []

for i in range(len(lau)):
	tmp = f[f.Latitude == lau[i]].reset_index(drop=True)
	tmp['CMAQ_NO2'] = tmp['CMAQ_NO2']*1000
	tmp.dt = pd.to_datetime(tmp.dt)
	#tmp.index = pd.to_datetime(tmp.dt)
	tmp.plot(x='dt',y='CMAQ_NO2',ax=ax[i+1],label='CMAQ',c='navy')
	tmp.plot(x='dt',y='Sample Measurement',ax=ax[i+1],label='EPA',c='royalblue',marker='o',linestyle='None',markersize=3)
	ax[i+1].set_title(str(i+1))
	ax[i+1].set_ylim([0,50])
	x,y = np.array(tmp['Sample Measurement']),np.array(tmp.CMAQ_NO2)
	stats.append([pearsonr(x,y)[0],np.mean(x-y),np.sum(x-y)/np.sum(x)])


pd.DataFrame(stats).to_csv('no2_perf.csv')
plt.tight_layout()
#ax[0].xaxis.set_major_formatter(date_form)
#fig.autofmt_xdate(rotation=90)
plt.savefig('NO2_prelim.pdf')
plt.show()

ll='/projects/b1045/wrf-cmaq/output/Chicago_LADCO/mcip/PXLSM/ChicagoLADCO_d03/latlon_ChicagoLADCO_d03.nc' 
lon,lat = np.array(Dataset(ll)['lon']),np.array(Dataset(ll)['lat'])

fig, ax = plt.subplots(subplot_kw={'projection': crs.PlateCarree()})
tiler = Stamen(style='terrain-background')
ax.add_image(tiler, 7)
ax.coastlines('10m')
ax.add_feature(cartopy.cfeature.LAND,facecolor='None',edgecolor='k')
states_provinces = cartopy.cfeature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='50m',facecolor='none')
ax.add_feature(states_provinces, edgecolor='gray')
[ax.text(x=lonu[i],y=lau[i],s=str(i+1)) for i in range(len(lonu))]

ax.set_xlim([lon.min(),lon.max()])
ax.set_ylim([lat.min(),lat.max()])
ax.set_title('NO2 Station Locations over 1.3 km Domain')
#plt.savefig('/projects/b1045/montgomery/paper1/counties_with_95_population.pdf')
plt.savefig('NO2_prelim_stnlocation.png',dpi=300)
plt.show()


################################################################################################################################################

f = pd.read_csv('O3_d03_2021_withCMAQ.csv',index_col=0)
f = f[f.x != 0].reset_index(drop=True)
f = f[f.x != 287].reset_index(drop=True)
f = f[f.y != 0].reset_index(drop=True)
f = f[f.y != 314].reset_index(drop=True)

f = f[~pd.isna(f['Sample Measurement'])].reset_index(drop=True)
x =  f['Sample Measurement']
y =  f['CMAQ_O3']

lau,lonu = f.Latitude.unique(),f.Longitude.unique()

fig,ax = plt.subplots(6,7,figsize=(16,9),sharex=True)
ax = ax.ravel()
plt.xticks(rotation=90)
stats = []

for i in range(len(ax)):
	tmp = f[f.Latitude == lau[i]].reset_index(drop=True)
	tmp['Sample Measurement'] = tmp['Sample Measurement']*1000
	tmp['CMAQ_O3'] = tmp['CMAQ_O3']*1000
	tmp.dt = pd.to_datetime(tmp.dt)
	#tmp.index = pd.to_datetime(tmp.dt)
	tmp.plot(x='dt',y='CMAQ_O3',ax=ax[i],c='navy')
	tmp.plot(x='dt',y='Sample Measurement',ax=ax[i],c='royalblue',marker='o',linestyle='None',markersize=3)
	ax[i].set_title(str(i))
	ax[i].set_ylim([0,100])
	x,y = np.array(tmp['Sample Measurement']),np.array(tmp.CMAQ_O3)
	stats.append([pearsonr(x,y)[0],np.mean(x-y),np.sum(x-y)/np.sum(x)])
	ax[i].get_legend().remove()

plt.tight_layout()
#ax[0].xaxis.set_major_formatter(date_form)
#fig.autofmt_xdate(rotation=90)
plt.savefig('o3_prelim_3.pdf')
plt.show()

ll='/projects/b1045/wrf-cmaq/output/Chicago_LADCO/mcip/PXLSM/ChicagoLADCO_d03/latlon_ChicagoLADCO_d03.nc' 
lon,lat = np.array(Dataset(ll)['lon']),np.array(Dataset(ll)['lat'])

fig, ax = plt.subplots(subplot_kw={'projection': crs.PlateCarree()})
tiler = Stamen(style='terrain-background')
ax.add_image(tiler, 9)
ax.coastlines('10m')
ax.add_feature(cfeature.LAND,facecolor='None',edgecolor='k')
states_provinces = cfeature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='50m',facecolor='none')
ax.add_feature(states_provinces, edgecolor='gray')
din.plot('NO2',ax=ax)
[ax.text(x=lonu[i],y=lau[i],s=str(i+1)) for i in range(len(lonu))]

ax.set_xlim([lon.min(),lon.max()])
ax.set_ylim([lat.min(),lat.max()])
ax.set_title('O3 Station Locations over 1.3 km Domain')
#plt.savefig('/projects/b1045/montgomery/paper1/counties_with_95_population.pdf')
plt.savefig('o3_prelim_stnlocation.png',dpi=350)
plt.show()



################################################################################################################################
# hourly plots -- 
# day vs night

season = 'Aug2021' #Aug2021

if season == 'Aug2021':
	f = pd.read_csv('/projects/b1045/montgomery/CMAQcheck/NO2_d03_2021_8_EPA_CMAQ_Combine.csv',index_col=0)
	ken = f[f['Site Num'] == 219]
	ec = pd.read_csv('/projects/b1045/montgomery/microsoft/cleaned_data_withCMAQ_0802021.csv')
	ec = ec[ec.DeviceId == 2136].reset_index(drop=True)
	ec['date'] = pd.to_datetime(pd.date_range("2021-08-01","2021-08-31 23:00",freq="H"))


if season == 'Feb2022':
	f = pd.read_csv('/projects/b1045/montgomery/CMAQcheck/NO2_d03_2022_2_EPA_CMAQ_Combine.csv',index_col=0)
	ken = f[f['Site Num'] == 219]
	ec = pd.read_csv('/projects/b1045/montgomery/microsoft/cleaned_data_wint.csv'); ec['date'] = pd.to_datetime(ec['Unnamed: 0'])
	ec = ec[ec.DeviceId == 2136].reset_index(drop=True)


ken['date'] = pd.to_datetime(ken["dt"]).reset_index(drop=True)

tmp1 = ec[['date','CalibratedNO2']]
tmp1.columns = ['dt','ec_no2']
tmp2 = ken[['dt','CMAQ','Sample Measurement']].reset_index(drop=True)
tmp2.dt = pd.to_datetime(tmp2.dt)
tmp2.columns = ['dt','cmaq_no2','epa_no2']

from dateutil import tz
from_zone = tz.gettz('UTC')
to_zone = tz.gettz('America/Chicago')
no2 = pd.merge(tmp1,tmp2,on='dt')
no2['h'] = [no2.dt[i].hour for i in range(len(no2))]
no2.index = no2.dt
no2 = no2.tz_localize(from_zone).tz_convert(to_zone)
no2 = no2.rename_axis('dt0').reset_index()
no2['h0'] = [no2.dt0[i].hour for i in range(len(no2))]

di = no2.groupby('h0').mean().reset_index()
dimin = no2.groupby('h0').std().reset_index()
#dimax = no2.groupby('h').max().reset_index()


################################################################
vs = ['ec_no2','cmaq_no2','epa_no2']
label=['Eclipse','WRF-CMAQ','EPA']
marker=['s','^','o']
off = [0,0.25,0.5]
colors = ["#FCAC60","#4E63AB","#4d934d"]#"#71846E"]#"#A80C44"]
colors = ["#2974bc","#4d934d","#FCAC60"]
#colors = ["#2c75ff","#FCAC60","#4d934d"]
sizes = [3.5,5,3.5]

fig,ax=plt.subplots(figsize=(4.5,4))

for i in range(len(vs)):
	ax.errorbar(di['h0']+off[i],di[vs[i]],yerr=dimin[vs[i]]/2,linestyle="None",linewidth=.9,marker=marker[i],markersize=sizes[i],label=label[i],c=colors[i])

ax.set_xlim([-0.2,24])
ax.set_xticks(np.arange(0,24)[::2])
#plt.legend()
plt.legend(bbox_to_anchor=(.95, 1.1),  borderaxespad=0, ncol = 3,fontsize=10)

ax.set_ylim(0,40)
ax.set_xlabel('Hour (CT)')
ax.set_ylabel(r'NO$_2$ (ppb)')
#ax.set_title('(a) August 2021')
ax.grid(alpha=0.2)
plt.tight_layout()
plt.savefig('%s_hourlyEPACMAQEclipse.png'%(season),dpi=300)
plt.show()


# repeat but normalized
################################################################

def normalize(di,v):
	di[v+'norm'] = (di[v] - di[v].min())/(di[v].max() - di[v].min())
	return di

di = normalize(di,'ec_no2')
di = normalize(di,'cmaq_no2')
di = normalize(di,'epa_no2')

no2 = normalize(no2,'ec_no2')
no2 = normalize(no2,'cmaq_no2')
no2 = normalize(no2,'epa_no2')

dimin = no2.groupby('h0').std().reset_index()

vs = ['ec_no2norm','cmaq_no2norm','epa_no2norm']


fig,ax=plt.subplots(figsize=(4.5,4))

for i in range(len(vs)):
	ax.errorbar(di['h0']+off[i],di[vs[i]],yerr=dimin[vs[i]]/2,linestyle="None",linewidth=.9,marker=marker[i],markersize=sizes[i],label=label[i],c=colors[i])

ax.set_xlim([-0.2,24])
ax.set_xticks(np.arange(0,24)[::2])
#plt.legend()

ax.set_ylim(-0.2,1.2)
ax.set_xlabel('Hour (CT)')
ax.set_ylabel(r'NO$_2$ (normalized)')
ax.set_title('(a) Aug 2021')
ax.grid()
plt.tight_layout()
plt.savefig('%s_hourlyEPACMAQEclipse_normalized.png'%(season),dpi=300)
plt.show()



stats_normalized(no2['epa_no2'],no2['ec_no2'])
rushhour = np.arange(6,10)
after = np.arange(4+12,6+12+1)
day = np.arange(6,18+1)
night = np.array(list(np.arange(19,24)) + list(np.arange(1,6)))


rush = no2[no2.h0.isin(rushhour)].reset_index(drop=True)
after = no2[no2.h0.isin(after)].reset_index(drop=True)
day = no2[no2.h0.isin(day)].reset_index(drop=True)
night = no2[no2.h0.isin(night)].reset_index(drop=True)
out = pd.DataFrame([stats_normalized(no2['epa_no2'],no2['cmaq_no2']),
stats_normalized(rush['epa_no2'],rush['cmaq_no2']),
stats_normalized(after['epa_no2'],after['cmaq_no2']),
stats_normalized(day['epa_no2'],day['cmaq_no2']),
stats_normalized(night['epa_no2'],night['cmaq_no2']),
stats_normalized(no2['epa_no2'],no2['ec_no2']),
stats_normalized(rush['epa_no2'],rush['ec_no2']),
stats_normalized(after['epa_no2'],after['ec_no2']),
stats_normalized(day['epa_no2'],day['ec_no2']),
stats_normalized(night['epa_no2'],night['ec_no2']),
])
out.columns = ['mu_d','mu_p','mb','nmb','rmse','nrmse','r']
out['data'] = ['CMAQ','CMAQ_morning','CMAQ_afternoon','CMAQ_day','CMAQ_night','Eclipse','Eclipse_morning','Eclipse_afternoon','Eclipse_day','Eclipse_night']
out.to_csv('Aug2021_hourlyEPACMAQEclipse.csv')




di = di.set_index('h0')
d2 = di.drop('h',axis=1)
d2.plot(marker='*',linestyle="None")

streets = gpd.read_file('/projects/b1045/montgomery/shapefiles/geo_export_28abe032-5cb4-4d93-9c9e-cf51d26e8169.shp')
# make meters EPS
fin = pd.read_csv('Average_Daily_Traffic_Counts.csv')
fin['geometry'] = [Point(fin.Longitude[i],fin.Latitude[i]) for i in range(len(fin))]
fin = gpd.GeoDataFrame(fin,crs = streets.crs)
fin = fin.to_crs(class1.crs)
fin['dist'] = fin.distance(single_track_line)



#http://apps.dot.illinois.gov/gist2/
from shapely.ops import unary_union, cascaded_union
import geopandas as gpd
fin = gpd.read_file('/projects/b1045/montgomery/microsoft/T2HWY2021.shp')
clip = gpd.read_file('/projects/b1045/montgomery/microsoft/all_data_grid_eclipse_rf.geojson')

def crop_return_shp(fin,target_crs,gdf_in):
   # have to: import geopandas as gpd
   # Warning -- will not work if you're cropping points?
   # fin = shapefile to crop to (a geopandas dataframe)
   # gdf_in = geopandas dataframe with your data 
   # target_crs = the crs you want to crop to (recommended to be the fin.crs or gdf_in.crs)
   # gdf_out = data cropped to shapefile
   # gdf_out_union = outside bounds of the gdf you just cropped
for i in range(1):
   fin = fin.to_crs(target_crs)
   gdf_in = gdf_in.to_crs(target_crs)
   fin_union = gpd.GeoSeries(unary_union(fin.geometry))
   gdf_out = gpd.clip(gdf_in,fin_union).reset_index(drop=True)
   # removing anything that was clipped weird
   gdf_out['dt'] = [(gdf_out.geometry[i].type != 'Point')  & (gdf_out.geometry[i].type != 'MultiLineString') & (gdf_out.geometry[i].type != 'LineString') for i in range(len(fin))]
   #gdf_out=gdf_out[fin.dt].reset_index(drop=True)
   gdf_out = gdf_out.reset_index(drop=True)
   fin_union = gpd.GeoDataFrame({'geometry':fin_union},crs=target_crs)
   #fin = fin[fin.geometry.type == 'Polygon'].reset_index(drop=True)
   return gdf_out,fin_union

fin,target_crs,gdf_in = clip,fin.crs,fin
hwy_clip = gdf_out
clip  = clip.to_crs(fin.crs)
hwy_clip,hwy_union = crop_return_shp(clip,fin.crs,fin)
hwy_clip.to_file('T2HWY2021_clip.shp')

streets = gpd.read_file('/projects/b1045/montgomery/shapefiles/geo_export_28abe032-5cb4-4d93-9c9e-cf51d26e8169.shp')
class1 = streets[streets['class'] == '1'].reset_index(drop=True)
class1 = class1.to_crs(hwy_clip.crs)

ok = gpd.sjoin(class1,hwy_clip)
ok.to_file('T2HWY2021_clip_merge.shp')

ok = gpd.read_file('/projects/b1045/montgomery/microsoft/shapefile/T2HWY2021_clip_merge.shp')

ok.plot('AADT',legend=True)
plt.show()


#----------------------------------------------------------
import cartopy.feature as cfeature
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy

central_lat = 37.5
central_lon = -96
extent = [-120, -70, 23.5, 50.5]
central_lon = np.mean(extent[:2])
central_lat = np.mean(extent[2:])

plt.figure(figsize=(3, 1.5))
ax = plt.axes(projection=ccrs.AlbersEqualArea(central_lon, central_lat))
ax.set_extent(extent)

ax.add_feature(cartopy.feature.OCEAN)
ax.add_feature(cartopy.feature.LAND, edgecolor='black',linewidth=0.5)
ax.add_feature(cartopy.feature.LAKES, edgecolor='black',linewidth=0.5)
ax.add_feature(cartopy.feature.RIVERS,alpha=0.5)
states_provinces = cfeature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='50m',facecolor='none',linewidth=0.5)
ax.add_feature(states_provinces, edgecolor='gray',linewidth=0.5,alpha=0.5)
ax.add_feature(cfeature.BORDERS.with_scale('50m'),linewidth=0.5,alpha=0.5)
ax.gridlines(linewidth=0.5,alpha=0.5)
#ax.scatter(-87.623177,41.881832,transform=ccrs.AlbersEqualArea(central_lon, central_lat))
plt.savefig('USA.png',dpi=300)


