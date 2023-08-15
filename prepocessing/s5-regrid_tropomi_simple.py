#regrid_tropomi_simple.py
import xesmf as xe
import pandas as pd
import xarray
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

def bilinear_interp(floncrop,flatcrop,mlon,mlat,data):
			# floncrop,flatcrop -- target grid
			# mlon,mlat -- input grid
			# data -- data to regrid
			# regrid data so that it can be added to metdata
			grid_out = {"lon": floncrop, "lat": flatcrop}
			grid_in = {"lon": mlon, "lat": mlat}
			regridder = xe.Regridder(grid_in, grid_out, "bilinear")
			gridout = regridder(data)
			return gridout

def find_index(stn_lon, stn_lat, wrf_lon, wrf_lat):
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
		x, y = np.where(c == np.min(c))
		#add indices of nearest wrf point station
		xx.append(x)
		yy.append(y)
	#
	xx=[xx[i][0] for i in range(len(xx))];yy=[yy[i][0] for i in range(len(yy))]
	#return indices list
	return xx, yy

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


# Step 1 do the regridding

fnames = list(pd.read_csv('f.csv')['0']) # names of files
fnames = np.array(pd.read_csv('f08.txt')).ravel() # names of files
#fnames = list(pd.read_csv('satnames_022022.csv')['0']) # names of files
flatcrop= np.array(pd.read_csv('cropped_chi_lat_ESPG4326.csv',index_col = 0)) # cmaq grid
floncrop= np.array(pd.read_csv('cropped_chi_lon_ESPG4326.csv',index_col = 0)) # cmaq grid

# regrid data so that it can be added to metdata
prec = []
no2_bi = []
lat = []
lon = []
prec_bi = []

for i in range(len(fnames)):
	try:
		#f = Dataset('./tropomi/'+fnames[i]+'.nc')['PRODUCT']
		f = Dataset(fnames[i])['PRODUCT']
		#f = Dataset('./tropomi/'+fnames[i])['PRODUCT']
	except:
		continue
	latmask = f['latitude'][0]
	lonmask = f['longitude'][0]
	precmask = f['qa_value'][0] # interpolate this so we preserve where missing data is
	data = f['nitrogendioxide_tropospheric_column'][0]
	# check time, remove if not near overpass time
	# make masks
	#mask = (latmask < 42) & (latmask > 41.5)  & (lonmask > -88) & (lonmask < -87.3)#& (precmask >= 0.75)
	maskloc = (np.where((latmask < 42.5) & (latmask > 41.5) & (lonmask > -88) & (lonmask < -87.3)))
	min_tuple = maskloc[0].min(),maskloc[1].min()
	max_tuple = maskloc[0].max(),maskloc[1].max()
	#mask to chi if overpass time is right
	t = int(f['time_utc'][0][min_tuple[0]].split('T')[1][:2])
	print(t)
	if t <= 19 & t >= 16:	
		data = data[min_tuple[0]:max_tuple[0],min_tuple[1]:max_tuple[1]]
		latmask = latmask[min_tuple[0]:max_tuple[0],min_tuple[1]:max_tuple[1]]
		lonmask = lonmask[min_tuple[0]:max_tuple[0],min_tuple[1]:max_tuple[1]]
		precmask = precmask[min_tuple[0]:max_tuple[0],min_tuple[1]:max_tuple[1]]
		# remove nans by filling with average value
		avg = np.nanmean(data)
		data = data.data
		data[data > 100000] = avg
		precmask = precmask.data
		precmask[precmask > 100000] = 0
		#data[precmask < 0.75] = avg
 		#prec.append(f['nitrogendioxide_tropospheric_column_precision_kernel'][0][mask])
		#interp to chi
		prec_check = bilinear_interp(floncrop,flatcrop,lonmask,latmask,precmask)
		if np.sum(prec_check < 0.75) < 600:
			no2_bi.append(bilinear_interp(floncrop,flatcrop,lonmask,latmask,data))
			prec_bi.append(prec_check)
			print('done with '+str(i))
		else:
			no2_bi.append(np.zeros(floncrop.shape)*np.nan)
			prec_bi.append(prec_check)
			print('done with '+str(i))
	# where doesnt this work
	else: print(i,t)

# For aug -- remove outlier
#no2_bi = np.array([no2_bi[i] for i in range(0,7)]+[no2_bi[i] for i in range(8,len(no2_bi))])
#prec_bi = np.array([prec_bi[i] for i in range(0,7)]+[prec_bi[i] for i in range(8,len(prec_bi))])

no2_bi_all = np.ma.masked_array(np.array(no2_bi),mask = np.array(prec_bi) < 0.75)
pixel_counts = np.sum(np.array(prec_bi) >= 0.75,axis=0)
#output pixel counts
#pd.DataFrame(pixel_counts).to_csv('pixel_counts_202202.csv')

labels = pd.read_csv('u_lon_lat_labels_for_grid.csv')
days = np.array([[i]*len(labels) for i in range(31)]).ravel()
out = pd.DataFrame([list(labels['u'])*31,days,no2_bi_all.filled(np.nan).ravel()]).T
out.to_csv('trop_daily_08.csv')

fig,ax = plt.subplots(6,6,figsize=(12,12))
ax=ax.ravel()

for i in range(len(no2_bi_all)):
	ax[i].pcolormesh(floncrop,flatcrop,no2_bi_all[i]*10**5,vmin=0,vmax=20,shading='auto')
	ax[i].set_title(fnames[i].split('T')[0].split('_')[-1]+'-'+fnames[i].split('_')[-1][-8:])

plt.tight_layout()
plt.show()





fig,ax = plt.subplots(8,5,figsize=(12,12))
ax=ax.ravel()

for i in range(len(no2_bi_all)):
	ax[i].pcolormesh(floncrop,flatcrop,no2_bi[i]*10**5,vmin=0,vmax=20,shading='auto')

plt.tight_layout()
plt.show()

fig,ax = plt.subplots(1,2,figsize=(12,6))
ax[0].pcolormesh(floncrop,flatcrop,np.nanmean(no2_bi_all*10**5,axis=0),vmin=5,vmax=25,shading='auto')
ax[1].pcolormesh(floncrop,flatcrop,np.nanmean(np.array(no2_bi)*10**5,axis=0),vmin=5,vmax=25,shading='auto')
plt.tight_layout()
plt.show()

#output 
#pd.DataFrame(np.nanmean(no2_bi_all*10**5,axis=0)).to_csv('02_2022_tropomi_NO2.csv')

# Step 2
# check perfromance of regridding

# pull in epa
epa = pd.read_csv('hourly_42602_2022.csv')
epa = epa[(epa['State Name'] == 'Illinois') & (epa['County Name'] == 'Cook')].reset_index(drop=True)
epa['dt'] = [epa['Date GMT'][i] + ' ' + epa['Time GMT'][i] for i in range(len(epa))]
epa.dt = pd.to_datetime(epa.dt)
#epa = epa[(epa.dt > pd.to_datetime('2021-08-01 00:00:00')) & (epa.dt < pd.to_datetime('2021-08-31 23:00:00'))]
epa = epa[(epa.dt > pd.to_datetime('2022-02-01 00:00:00')) & (epa.dt < pd.to_datetime('2022-02-28 23:00:00'))]
epa = epa.groupby('Site Num').mean().reset_index()

indx,indy = find_index(epa.Longitude.unique(),epa.Latitude.unique(),floncrop,flatcrop)
#indx0,indy0 = indx[1:],indy[1:]
epa_lon, epa_lat = epa.Longitude.unique(),epa.Latitude.unique()

epa_monthly = np.array(epa['Sample Measurement'])
#no2_bi = np.array(no2_bi)*10**5
#no2_monthly_trop = np.nanmean(no2_bi,axis=0)*10**5
trop_monthly = np.array(np.nanmean(no2_bi,axis=0)[indx,indy])
trop_monthly = trop_monthly*10**5
epa_stat = stats_normalized(epa_monthly,trop_monthly)

# pull in eclipse
fin = pd.read_csv('ec_on_cmaq.csv',index_col = 0)
#fin = gpd.read_file('all_data_grid_eclipse_rf.geojson')
avg_obs = fin.groupby('u').mean().reset_index()

# check against eclipse and eclipse 
ec = list(avg_obs[~np.isnan(avg_obs.ec_NO2)].ec_NO2) 
eclat,eclon = list(avg_obs[~np.isnan(avg_obs.ec_NO2)].lat),list(avg_obs[~np.isnan(avg_obs.ec_NO2)].lon) 
indx1,indy1 = find_index(eclon,eclat,floncrop,flatcrop)
#for feb == pull ec
fin = gpd.read_file('all_data_grid_eclipse_rf.geojson')
latlon_idx = [np.argmin(np.maximum(np.abs(fin.x-eclon[i]),np.abs(fin.y-eclat[i]))) for i in range(len(eclon))]

ec_stat = stats_normalized(fin['rf_month_february'][latlon_idx],np.array(no2_monthly_trop[indx1,indy1]))

# check against eclipse and eclipse + EPA
ecepa = list(epa_monthly)+list(ec)
indx,indy = find_index(epa.Longitude.unique(),epa.Latitude.unique(),floncrop,flatcrop)
indx,indy = indx+indx1,indy+indy1
ec_epa_stat = stats_normalized(np.array(ecepa),np.array(no2_monthly_trop[indx,indy]))

ec = np.array(ec)

# do normalized stats
nepa_monthly = (epa_monthly-np.min(epa_monthly))/(np.max(epa_monthly) - np.min(epa_monthly))
ntrop_monthly =(no2_monthly_trop-np.min(no2_monthly_trop))/(np.max(no2_monthly_trop) - np.min(no2_monthly_trop))
necepa = (ecepa-np.min(ecepa))/(np.max(ecepa) - np.min(ecepa))
nec =  (ec-np.min(ec))/(np.max(ec) - np.min(ec))

epa_stat = stats_normalized(nepa_monthly,ntrop_monthly[indx0,indy0])
ec_stat = stats_normalized(np.array(nec),np.array(ntrop_monthly[indx1,indy1]))
ec_epa_stat = stats_normalized(np.array(necepa),np.array(ntrop_monthly[indx,indy]))


stats = pd.DataFrame([epa_stat,ec_stat,ec_epa_stat])
stats.columns = ['n_mu_d','n_mu_p','mb','nmb','rmse','nrmse','r']
stats.to_csv('tropomi_vs_epa_ec_08.csv')


#
tmp_shp = gpd.read_file('rectangle_chicago_1km_grid.shp')
tmp_shp['trop_02'] = np.array(pd.read_csv('02_2022_tropomi_NO2.csv',index_col = 0)).ravel()
tmp_shp['trop_08'] = np.array(pd.read_csv('08_2021_tropomi_NO2.csv',index_col = 0)).ravel()
tmp_shp['trop_02_ct'] = np.array(pd.read_csv('pixel_counts_202108.csv',index_col = 0)).ravel()
tmp_shp['trop_08_ct'] = np.array(pd.read_csv('pixel_counts_202202.csv',index_col = 0)).ravel()
tmp_shp['cmaq_08'] = np.array(pd.read_csv('cropped_no2_monthly_ESPG4326.csv',index_col = 0)).ravel()
tmp_shp['cmaq_02'] = np.array(pd.read_csv('cropped_no2_monthly_feb2022_ESPG4326.csv',index_col = 0)).ravel()

merged = gpd.GeoDataFrame(pd.merge(tmp_shp,fin,on='u'),crs = fin.crs)
merged['geometry'] = merged.geometry_y

import matplotlib
fig,ax=plt.subplots(2,3)
ax=ax.ravel()
my_cmap = copy.copy(matplotlib.cm.get_cmap('gray')) # copy the default cmap
my_cmap.set_bad((0,0,0))
titles=['CMAQ','Eclipse','TropOMI']
cols = ['cmaq_08','rf_month_august','trop_08','cmaq_02','rf_month_february','trop_02']
vmins = [5,5,5,5,5,5]
vmaxs = [20,20,20,20,20,20]
for i in range(len(ax)):
	merged.plot(cols[i],ax=ax[i],vmin=1,vmax=20,norm=matplotlib.colors.LogNorm(vmin=5,vmax=25),
           legend=True)
	ax[i].axis('off')
	if i < 3: ax[i].set_title(titles[i])
	if i == 0: ax[i].set_ylabel('Aug. 2021')
	if i == 3: ax[i].set_ylabel('Feb. 2022')
	print(merged[cols[i]].mean())

plt.tight_layout()
#plt.savefig('datasets.png',dpi=300)
plt.show()

####

fig,axs = plt.subplots(2,3,figsize=(12,9))
axs = axs.ravel()

names = ['cmaq','tropomi','eclipse']*2
titles = ['CMAQ\n(a)','TropOMI\n(b)','Eclipse\n(c)','(d)','(e)','(f)']
vmins,vmaxs = [5,5,5]*2,[15,15,15,20,20,20]
units = ['ppb',r'molecules/cm$^2$*10$^5$','ppb']*2

for i in range(0,3):
    ax = axs[i]
    dgs.plot(names[i],ax=ax,facecolor='k')
    c1=dgs.plot(names[i],ax=ax,vmin = vmins[i],vmax = vmaxs[i])
    class1.plot(ax=ax,edgecolor='k',alpha = 0.3)
    ax.text(-0.05,-0.01,s=r'$\mu$ = %.1f %s'%(dgs[names[i]].mean(),units[i]),transform=ax.transAxes)
    ax.set_title(titles[i])
    ax.axis('off')

for i in range(3,6):
    ax = axs[i]
    dgw.plot(ax=ax,facecolor='k')
    c2=dgw.plot(names[i],ax=ax,vmin = vmins[i],vmax = vmaxs[i])
    class1.plot(ax=ax,edgecolor='k',alpha = 0.3)
    ax.set_title(titles[i],fontsize=12)
    ax.text(-0.05,-0.01,s=r'$\mu$ = %.1f %s'%(dgw[names[i]].mean(),units[i]),transform=ax.transAxes)
    ax.axis('off')

# add colorbar
cax = fig.add_axes([0.9, 0.55, 0.02, 0.3]) #locaation in plot -x ,starting point from bottom-y ,width,length
sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=vmins[0], vmax=vmaxs[0]))
sm._A = []
fig.colorbar(sm, cax=cax,label = r'ppb, molecules/cm$^2$ * 10$^5$')

# add colorbar
cax = fig.add_axes([0.9, 0.15, 0.02, 0.3]) #locaation in plot -x ,starting point from bottom-y ,width,length
sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=vmins[0], vmax=vmaxs[0]))
sm._A = []
fig.colorbar(sm, cax=cax,label = r'ppb, molecules/cm$^2$ * 10$^5$')

# add label
axs[0].text(-0.1,0.5,s='Aug. 2021',transform=axs[0].transAxes,rotation=90,fontsize = 12)
axs[3].text(-0.1,0.5,s='Feb. 2022',transform=axs[3].transAxes,rotation=90,fontsize = 12)

plt.show()


####

titles = ['winddir_re','airt_rean','rh_rean','ws_rean','ndvi','arterials','arterials','industry_a','res_ar','com_ar','manu_ar','open_ar']
units = ['degree','K','%','m/s',' ','counts','counts','fraction','fraction','fraction','fraction','fraction','fraction']
vmins = [150,301.5,43,3.5]+[0]*11
vmaxs = [180,302.5,53,4.5,1,60,60]+[1.0]*6

fig,axs = plt.subplots(4,3,figsize=(8,10))
axs = axs.ravel()

for i in range(len(axs)):
	ax = axs[i]
	if i == 5: ax.axis('off')
	else:
		axi = dg.plot(titles[i],ax=ax,facecolor='k',vmin=vmins[i],vmax=vmaxs[i],legend=True,legend_kwds={"shrink":.5,"label":units[i]})
		c = class1.plot(ax=ax,edgecolor='k',alpha = 0.3)
		ax.text(-0.05,-0.01,s=r'$\mu$ = %.1f'%(dg[titles[i]].mean()),transform=ax.transAxes)
		ax.set_title(titles[i])
		ax.axis('off')

plt.tight_layout()
plt.savefig('chicago_chars.png',dpi=300)
plt.show()

# add colorbar
cax = fig.add_axes([0.9, 0.15, 0.02, 0.3]) #locaation in plot -x ,starting point from bottom-y ,width,length
sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=vmins[0], vmax=vmaxs[0]))
sm._A = []
fig.colorbar(sm, cax=cax,label = r'ppb, molecules/cm$^2$ * 10$^5$')

# add label
axs[0].text(-0.1,0.5,s='Aug. 2021',transform=axs[0].transAxes,rotation=90,fontsize = 12)
axs[3].text(-0.1,0.5,s='Feb. 2022',transform=axs[3].transAxes,rotation=90,fontsize = 12)

plt.show()
