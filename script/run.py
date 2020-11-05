#!/usr/bin/env python
# coding: utf-8

__author__ = 'Emergency Division WFP'
__project__ = 'FbF project Maguindanao, WFP Philippine CO'
__contact__ = 'michael.manalili@wfp.org', 'wfp.hq.gis@wfp.org'

#!/usr/bin/env python3
import requests, os, time, shutil
from time import strptime
import matplotlib.pyplot as plt
from datetime import datetime, timedelta, date
import pandas as pd
import folium
import geopandas as gpd
import descartes
import schedule
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart #python3
from email.mime.base import MIMEBase #python3
from email import encoders #python3
from sentinelsat.sentinel import SentinelAPI, read_geojson, geojson_to_wkt
import ee
import eeconvert
import rasterio
from rasterio.plot import show
from rasterstats import zonal_stats
import cdsapi
from geopandas import GeoSeries, GeoDataFrame
from sqlalchemy import create_engine
import psycopg2 
import io
from lxml import html
import wget
import cdsapi
from configparser import ConfigParser
import sentinelsat
from sqlalchemy import create_engine
import psycopg2 
import pandas.io.sql as psql
from io import BytesIO as StringIO
from osgeo import ogr
from ftplib import FTP
import pycountry
import uuid
#import pysftp
#from io import BytesIO as StringIO

try:
    ee.Initialize()
    print ('Earth Engine package initialized successfully..')
except ee.EEException as e:
    print ('Earth Engine package failed to initialize!')
except:
    print ('Unexpected error:', sys.exc_info()[0])
    raise

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

pd.set_option('mode.chained_assignment', None)
#pd.set_option('mode.chained_assignment', 'raise')

execTime = time.time()

now = datetime.today()
d = str(now)
date = d[0:10]
date_cut = date.replace('-', '')
whitelist_data = 'extreme_precip_alert_' + date_cut
alert_data = 'red_alert_' + date_cut

script = os.path.dirname(os.path.realpath("__file__"))
root = os.path.dirname(script)

# Assumed this folders are pre-created
data_root = root + '/data/'
process_root = root + '/data/processing'
process_ecmwf = process_root + '/process_ecmwf/'

alert_csv_path = root + '/data/alert/'
map_view_path = root + '/data/view/'

alert_data_name = 'FbF_Targetting_PHL_' + date_cut + '.xlsx'
map_view_name = 'Map_view_' + date_cut + '.html'

config = ConfigParser()
config.read(os.path.join(root, "config.txt"))

esa_un = config.get('sentinelsat', 'ss_un')
esa_pw = config.get('sentinelsat', 'ss_pw')
esa_link = config.get('sentinelsat', 'ss_link')
api = SentinelAPI(esa_un, esa_pw, esa_link)
api

config.read(os.path.join(root, "config.txt"))
smtp_un = config.get('dbtrack','dbt_uname')
smtp_pw = config.get('dbtrack','dbt_pword')
print('config loaded.')

wfp_adm0 = ee.FeatureCollection("projects/unwfp/HQGIS/fbf/PHL_Maguindanao_adm4_pov")
data = eeconvert.fcToGdf(wfp_adm0)
l = list(data)
l.remove('UID')

def spatial_red(ans, reducer, scale):
    data = ans.reduceRegions(
            collection = wfp_aoi,
            reducer = reducer,
            crs='EPSG:3857',
            scale=scale,
            tileScale=1)
    return data   

def write2file(filename, header, gdf):
    data = gdf.to_csv(filename, columns=header, encoding='utf-8')
    return data

def export2drive_raster(raster, bbox ,prefix,scale):
    run = ee.batch.Export.image.toDrive(
         image=raster, 
         scale=scale,
         fileNamePrefix = prefix,
         region = ee.Feature(bbox).geometry().bounds().getInfo()['coordinates'], #can also be single country
         maxPixels = 1e13,
         fileFormat='GeoTiff',
         cloudOptimized=True).start()
    return run

#Format Options: GeoJSON, SHP, KML
def export2drive_table(ee_featCol, description):
    table = ee.batch.Export.table.toDrive(
         collection=ee_featCol, 
         description=description,
         fileFormat='CSV')#.start() 
    return table

def month(n):
    start = ee.Date('2014-01-01').advance(n, 'month')
    end = start.advance(1,'month')
    return ee.ImageCollection('JAXA/GPM_L3/GSMaP/v6/operational').select('hourlyPrecipRate').filterDate(start,end).mean().set('system:time_start', start.millis())

def get_wc():
    #x = ee.Image(collectionList.get(band_month))
    wc = ee.Image('users/wfphqgis/CLIM/wc20_30s_prec_'+s_bm)
    wc_vector_mean = spatial_red(wc,ee.Reducer.median(),900)
    df_wc = eeconvert.fcToGdf(wc_vector_mean)
    dff_wc = df_wc.drop(columns = l)
    dff_wc_r = dff_wc.rename(columns={'median':'WC_Mean_Prec'})
    return dff_wc_r

def get_sm2map():
    #x = ee.Image(collectionList.get(band_month))
    sm2 = ee.Image('users/wfphqgis/CLIM/clm_precipitation_sm2rain_m_1km_s00cm_20072018_v02_'+s_bm)
    sm2_vector_mean = spatial_red(sm2,ee.Reducer.median(),900)
    df_sm2 = eeconvert.fcToGdf(sm2_vector_mean)
    dff_sm2 = df_sm2.drop(columns = l)
    dff_sm2_r = dff_sm2.rename(columns={'median':'Mean_Prec'})
    return dff_sm2_r

def get_ecmwf():
    data_ecmwf = ee.Image('users/wfphqgis/CLIM/ecmwf_reanalysis_ymonmeanR_ymonmean')
    crs = data_ecmwf.projection().crs()
    scale = 9000
    tp = data_ecmwf.select('b' + s_bm)
    ecmwf = tp.resample('bilinear').reproject(crs = crs, scale = scale)
    wc_vector_mean = spatial_red(tp,ee.Reducer.median(),9000)
    df_wc = eeconvert.fcToGdf(wc_vector_mean)
    dff_wc = df_wc.drop(columns = l)
    dff_wc_r = dff_wc.rename(columns={'median':'ECMWF_Mean_Prec'})
    return dff_wc_r

def get_amsr2():
    #data = 'http://www.gdacs.org/flooddetection/DATA/AMSR2/AvgSignalTiffs/2019' 400K max
    data = 'http://www.gdacs.org/flooddetection/DATA/AMSR2/MagTiffs/2019' #20K max
    now = datetime.today()
    d = str(now)
    date = d[0:10]
    date_frmt = date.replace('-', '')

    page = requests.get(data)
    webpage = html.fromstring(page.content)

    amsr2 = webpage.xpath('//a/@href')
    amsr2_last = amsr2[-1]
    cur_year = d[0:4]
    base = 'http://www.gdacs.org'
    #prod = 'signal_4days_avg_4days_'
    prod = 'mag_signal_'
    hd5 = '_HD5'
    ext = '.tif'
    lnk =  base + amsr2_last + prod + date_frmt +hd5 +ext
    dl = wget.download(lnk,process_amsr2)
    return dl

def get_s1_info(i):
    s1 = api.query(final_df.envelope[i],
                   date=sensing_date,
                   platformname='Sentinel-1',
                   producttype='GRD',
                   area_relation='Intersects') #Intersects, IsWithin
                   #orbitdirection='ASCENDING')
    s1_json = api.to_geojson(s1)
    s1_gdf = api.to_geodataframe(s1)
    s1_gdf

    s1_gdf = s1_gdf[['ingestiondate']]#,'beginposition','ingestiondate','platformname','orbitnumber', 'producttype', 'endposition']]
    s1_gdf.sort_values('ingestiondate', ascending=True).head()
    s1_gdf

    dates = s1_gdf["ingestiondate"].dt.date.values
    map_date = dates + timedelta(days=11)
    ld = len(dates)
    
    print_path = open(root + '/data/alert/' + 'sentinel_sched_' + date_cut + '.txt','a')
    
    p = print("Sentinel Acquisition Schedule: " + "{0}: {1} \n{2} alert! \nThere are {3} Sentinel-1 overpass \nMappable on:\n{4}\n"
             .format(final_df['adm0_name'][i],final_df['adm2_name'][i],final_df['status'][i],ld,map_date),
             file=print_path)
    return p

def job():
    msg = MIMEMultipart()
    msg['From'] = 'adam.floods.wfp@gmail.com'
    msg['Subject'] = "Geographical Prioritization for FbF - WFP Philippines " + date_cut
    
    body = """<div style="color:rgb(28, 130, 196)"><font size="+2"><b><img id="myimage" src="https://mw1.google.com/crisisresponse/icons/un-ocha/activity_deployment_64px_icon_bluebox.png"></img>
    <br /> Geographical Prioritization for FbF - WFP Philippines</b> </font></div>
    <br/><b>Forecast run for %s. Please See attached summary table and map viewer.</b></b>
    <br />
    <br/><b>Note: </b></b>
    <br/>The information contained in this table is a guide only based on best available rainfall forecast information.
    <br />
    <br/>This is an automatically generated email, please do not reply.<br/>
    <br/>--------------------------------------------------------------<br/>
    <br/>Service provided by 
    <br/><i>UNITED NATIONS WORLD FOOD PROGRAMME
    <br/>WFP Emergency Division
    <br/>Contact: <a href="mailto:hq.gis@wfp.org">HQ Geospatial Support Unit</a></i></p>
    """ %(date_cut)#add 'a' in the modulo

    alert_fn = alert_csv_path + alert_data_name
    mapview_fn = map_view_path + map_view_name
    
    alerts = MIMEBase('application', "octet-stream")
    mapview = MIMEBase('application', "octet-stream")
    
    alerts.set_payload(open(alert_fn, "rb").read())
    mapview.set_payload(open(mapview_fn, "rb").read())
    
    encoders.encode_base64(alerts)
    encoders.encode_base64(mapview)
    
    alerts.add_header('Content-Disposition', 'attachment; filename=%s'%alert_data_name)
    mapview.add_header('Content-Disposition', 'attachment; filename=%s'%map_view_name)
    
    msg.attach(alerts)
    msg.attach(mapview)
    
    msg.attach(MIMEText(body, 'html'))
    
#     s = smtplib.SMTP('smtp.gmail.com', 587) #Anywhere
#     s.starttls() #Anywhere
    
    s = smtplib.SMTP_SSL('smtp.gmail.com', 465) # WFP domain
    
    s.ehlo()
    s.login(smtp_un, smtp_pw)
    sender = 'wfp.hq.dbtrack@gmail.com'
    
#     recipients = ['michaelandrew.manalili@gmail.com', 'michael.manalili@wfp.org' ,'mishael.argonza@wfp.org','damien.fontaine@wfp.org','martin.parreno@wfp.org',
#                   'jesse.mason@wfp.org','joan.odena@wfp.org','juanito.berja@wfp.org','isabelle.lacson@wfp.org','lara.prades@wfp.org','paris.kazis@wfp.org','abdel-lathif.younous@wfp.org']
    
    recipients = ['michaelandrew.manalili@gmail.com']
    
    if gbl_alert.empty == True:
        pass
        print('DataFrame is Empty, skipping email broadcast...')
        
    else:
        s.sendmail(sender, recipients, str(msg))
        s.quit()
        print('ALERT sent to subscribers!')

#dfo = GBL_dfo.filterBounds(wfp_adm0)
wfp_aoi = ee.FeatureCollection("projects/unwfp/HQGIS/fbf/PHL_Maguindanao_adm4_pov")
#wfp_aoi = GBL_adm2.filterBounds(dfo)
wfp_adm0 = ee.FeatureCollection("projects/unwfp/HQGIS/fbf/Maguindanao")
#dfodata = eeconvert.fcToGdf(GBL_dfo)


GSMaP = ee.ImageCollection('JAXA/GPM_L3/GSMaP/v6/operational')
precipitation = GSMaP.select('hourlyPrecipRate')

worldpop = ee.ImageCollection('WorldPop/GP/100m/pop').filter(ee.Filter.equals('year', 2020))
wpop2020 = worldpop.select('population').reduce(ee.Reducer.max()).clip(wfp_adm0)

ciesn2020 = ee.ImageCollection('CIESIN/GPWv4/unwpp-adjusted-population-count')
ciesn2020_pop = ciesn2020.select('population-count').reduce(ee.Reducer.max()).clip(wfp_adm0)

jrc_1m_depth = ee.Image('projects/unwfp/HQGIS/fbf/Maguindanao_FL_haz_bin')

hectares = ee.Image.pixelArea().divide(10000)

hectares_jrc = jrc_1m_depth.multiply(hectares)

pop_haz_area = wpop2020.multiply(jrc_1m_depth)
#pop_haz_area = ciesn2020_pop.multiply(jrc_1m_depth)

#Change
pxr = ee.Image('users/wfphqgis/CLIM/Global_IDF')
pxr_3days = pxr.select('b12') #band 12 for 72hours (3Days) band 14 is 120hours (5days) and 
pxr_3days_vector = spatial_red(pxr_3days,ee.Reducer.mode(),30000)

#Uncomment for GFS data (update datetime is not yet fixed)
#gfs_fc = ee.ImageCollection('NOAA/GFS0P25').select("total_precipitation_surface").filter(ee.Filter.eq('creation_time',ee.Date(0).update(2019,7,21,6,0,0)
#                                                                                                      .millis())).filter(ee.Filter.eq('forecast_hours',72)).sum()
#gfs_fc72h_res = gfs_fc.resample('bilinear').reproject(crs=gfs_fc.projection(), scale=10000).clip(wfp_adm0)

#GEE Datetime for GPM accumulated Precipitation
day0 = datetime.today()
hrs24 = day0 - timedelta(days=1)
hrs48 = day0 - timedelta(days=2)
hrs72 = day0 - timedelta(days=3)
day7 = day0 - timedelta(days=7)
day10 = day0 - timedelta(days=10)

#Temporal Reducer
gpm_precip_24h = precipitation.filter(ee.Filter.date(hrs24, day0))
gpm_precip_48h = precipitation.filter(ee.Filter.date(hrs48, day0))
gpm_precip_72h = precipitation.filter(ee.Filter.date(hrs72, day0))
gpm_precip_day7 = precipitation.filter(ee.Filter.date(day7, day0))
gpm_precip_day10 = precipitation.filter(ee.Filter.date(day10, day0))

#Changed to WorldClimV2
band_month = datetime.today().strftime("%m")
#band_month = strptime(run_date,'%b').tm_mon
s_bm = str(band_month)
dff_wc_rename = get_wc()
#dff_ecmwf_rename = get_ecmwf()
dff_sm2_rename = get_sm2map()
print('Climate Data loaded..')

#Spatial Reducer
gpm72h = gpm_precip_72h.sum().rename('72h').clip(wfp_aoi)
gpm10D = gpm_precip_day10.sum().rename('10D').clip(wfp_aoi)
gsmap_vector_mean = spatial_red(gpm72h,ee.Reducer.mode(),500)
gsmap_vector_mean2 = spatial_red(gpm10D,ee.Reducer.mode(),500)

wpop_vector = spatial_red(wpop2020,ee.Reducer.sum(),95)
jrc_vector = spatial_red(jrc_1m_depth,ee.Reducer.sum(),1000) #1km resolution all models
exp_vector = spatial_red(pop_haz_area,ee.Reducer.sum(),1000) #95 for WorldPop and 1000 for CIESN

#Change
ciesn_vector = spatial_red(ciesn2020_pop,ee.Reducer.sum(),1000)
#gfs_fc72h_vector = spatial_red(gfs_fc72h_res,ee.Reducer.mode(),10000)


print('Computing Global parameters. Please wait...')
df_precip_m = eeconvert.fcToGdf(gsmap_vector_mean)
df_precip_10 = eeconvert.fcToGdf(gsmap_vector_mean2)

df_pop = eeconvert.fcToGdf(wpop_vector)
df_jrc = eeconvert.fcToGdf(jrc_vector)
df_exp = eeconvert.fcToGdf(exp_vector)

#Change
df_pxr = eeconvert.fcToGdf(pxr_3days_vector)
df_pop1km = eeconvert.fcToGdf(ciesn_vector)
#df_gfs = eeconvert.fcToGdf(gfs_fc72h_vector)

#ECMWF datetime format for FTP (Temporary Solution while migrating to Windows)
forecast = datetime.today()
fc_fmt = forecast.strftime("%m%d")
days3 = forecast + timedelta(days=3)
days5 = forecast + timedelta(days=5)
days7 = forecast + timedelta(days=7)
days10 = forecast + timedelta(days=10)

days3F = days3.strftime("%m%d")
days5F = days5.strftime("%m%d")
days7F = days7.strftime("%m%d")
days10F = days10.strftime("%m%d")

ecmwf_3D_fc = fc_fmt + days3F + '00' + 'TP' + '.tiff'
ecmwf_5D_fc = fc_fmt + days5F + '00' + 'TP' + '.tiff'
ecmwf_7D_fc = fc_fmt + days7F + '00' + 'TP' + '.tiff'
ecmwf_10D_fc = fc_fmt + days10F + '00' + 'TP' + '.tiff'

config.read(os.path.join(root, "config.txt"))
ftp_user = config.get('gisftp', 'ftp_un')
ftp_pw = config.get('gisftp', 'ftp_pw')
ftp_gis =  config.get('gisftp', 'ftp_url')
ftp = FTP(ftp_gis)
ftp.login(user=ftp_user, passwd = ftp_pw)
save2local = ftp.cwd("/ECMWF_processed")

def grabFile(fn):
    filename = fn
    localfile = open(process_ecmwf + filename, 'wb')
    return ftp.retrbinary('RETR ' + filename, localfile.write, 1024)


grabFile(ecmwf_3D_fc)
grabFile(ecmwf_5D_fc)
grabFile(ecmwf_7D_fc)
grabFile(ecmwf_10D_fc)

ftp.quit()

dff_pop = df_pop.drop(columns = ['ADM1_PCODE', 'ADM2_PCODE', 'ADM3_PCODE', 'Barangay', 'Province', 'Region'])
dff_pop_rename = dff_pop.rename(columns={'sum':'Popn_VUL_2'})

dff_precip = df_precip_m.drop(columns = l)
dff_p_rename = dff_precip.rename(columns={'mode':'3Days_AccRF'})

dff_precip10 = df_precip_10.drop(columns = l)
dff_p_rename10 = dff_precip10.rename(columns={'mode':'10Days_AccRF'})

dff_jrc = df_jrc.drop(columns = l)
dff_jrc_rename = dff_jrc.rename(columns={'sum':'Flood_HAZ'})

dff_exp = df_exp.drop(columns = l)
dff_exp_rename = dff_exp.rename(columns={'sum':'Popn_at_Risk_FL'})

#Change PXR RFIDF
dff_pxr = df_pxr.drop(columns = l)
dff_pxr_rename = dff_pxr.rename(columns={'mode':'Max_3D_TP'})

dff_pop1km = df_pop1km.drop(columns = l)
dff_pop1km_rename = dff_pop1km.rename(columns={'sum':'Popn_VUL_1'})

#Merge all DF
pre_merged = dff_p_rename.merge(dff_sm2_rename, on='UID').merge(dff_pop_rename,on='UID').merge(dff_jrc_rename,on='UID').merge(dff_exp_rename,on='UID').merge(dff_pxr_rename,on='UID').merge(dff_pop1km_rename,on='UID').merge(dff_p_rename10,on='UID')#.merge(df_gfs_rename,on='adm2_code')
pre_merged['Max_3D_TP'] = pre_merged['Max_3D_TP'] * 1000
pre_merged = pre_merged.round(3)
pre_merged.tail(2)

pre_merged = gpd.GeoDataFrame(pre_merged, geometry=pre_merged['geometry'])

### Uncomment for NRT
ecmwf_F_3D = process_ecmwf + ecmwf_3D_fc
ecmwf_F_5D = process_ecmwf + ecmwf_5D_fc
ecmwf_F_7D = process_ecmwf + ecmwf_7D_fc
ecmwf_F_10D = process_ecmwf + ecmwf_10D_fc

ecmwf_F_3D_stats = zonal_stats(pre_merged, ecmwf_F_3D, prefix='3D_RF_forecast_',geojson_out=True)
ecmwf_F_5D_stats = zonal_stats(pre_merged, ecmwf_F_5D, prefix='5D_RF_forecast_',geojson_out=True)
ecmwf_F_7D_stats = zonal_stats(pre_merged, ecmwf_F_7D, prefix='7D_RF_forecast_',geojson_out=True)
ecmwf_F_10D_stats = zonal_stats(pre_merged, ecmwf_F_10D, prefix='10D_RF_forecast_',geojson_out=True)

ecmwf_F_3D_gdf = GeoDataFrame.from_features(ecmwf_F_3D_stats)
ecmwf_F_5D_gdf = GeoDataFrame.from_features(ecmwf_F_5D_stats)
ecmwf_F_7D_gdf = GeoDataFrame.from_features(ecmwf_F_7D_stats)
ecmwf_F_10D_gdf = GeoDataFrame.from_features(ecmwf_F_10D_stats)

#RENAME ME
ecmwf_F_3D_gdf['3D_RF_forecast'] = ecmwf_F_3D_gdf['3D_RF_forecast_mean'] * 1000 # converts TP values from meters to mm
ecmwf_F_5D_gdf['5D_RF_forecast'] = ecmwf_F_5D_gdf['5D_RF_forecast_mean'] * 1000 # converts TP values from meters to mm
ecmwf_F_7D_gdf['7D_RF_forecast'] = ecmwf_F_7D_gdf['7D_RF_forecast_mean'] * 1000 # converts TP values from meters to mm
ecmwf_F_10D_gdf['10D_RF_forecast'] = ecmwf_F_10D_gdf['10D_RF_forecast_mean'] * 1000 # converts TP values from meters to mm

dfNew = ecmwf_F_5D_gdf.merge(ecmwf_F_3D_gdf, left_index=True, right_index=True,how='outer', suffixes=('', '_y')).merge(ecmwf_F_7D_gdf, left_index=True, right_index=True,how='outer', suffixes=('', '_y')).merge(ecmwf_F_10D_gdf, left_index=True, right_index=True,how='outer', suffixes=('', '_y'))
dfNew.drop(list(dfNew.filter(regex='_y$')), axis=1, inplace=True)
to_drop = ['3D_RF_forecast_count','3D_RF_forecast_max','3D_RF_forecast_min','3D_RF_forecast_mean','5D_RF_forecast_count','5D_RF_forecast_max','5D_RF_forecast_min','5D_RF_forecast_mean','7D_RF_forecast_count','7D_RF_forecast_max','7D_RF_forecast_min','7D_RF_forecast_mean','10D_RF_forecast_count','10D_RF_forecast_max','10D_RF_forecast_min','10D_RF_forecast_mean']
gbl_alert = dfNew.drop(columns=to_drop).round(2)
gbl_alert['uuid'] = [uuid.uuid4() for _ in range(len(gbl_alert.index))]
print('Global Alerts created...')
#gbl_alert = gbl_alert.loc[(gbl_alert['3D_acc_precip'] >= 50)]

def percent_normal(x):
    if x <= 40:
        return "Way Below Normal"
    elif 41 <= x <= 80:
        return "Below Normal"
    elif 81 <= x <= 120:
        return "Near Normal"
    else:
        return "Above Normal"
    
def antecedent(x):
    if x <= 100:
        return "Light"
    elif 100 <= x <= 150:
        return "Moderate"
    elif 150 <= x <= 200:
        return "Intense"
    elif 200 <= x <= 250:
        return "Heavy"
    else:
        return "Torrential"
    
final_df = gbl_alert

#PAGASA Mean Rainfall for November
nov_monthly_rainfall_ave = 216

#See notes for values
final_df['running_rainfall_pn'] = ((final_df['10Days_AccRF'] + final_df['10D_RF_forecast'])/nov_monthly_rainfall_ave) * 100
final_df['NOV_%Normal_RF'] = final_df['running_rainfall_pn'].apply(percent_normal)

fbf_sug = "FbF_suggestion_for " + date_cut
final_df[fbf_sug] = ""

rf_cat = "RF_Forecast_cat_" + date_cut
final_df[rf_cat] = ""

# #See notes for values
final_df[rf_cat] = (final_df['3Days_AccRF'] +  final_df['10D_RF_forecast']).apply(antecedent)

# #See notes for values
# final_df.loc[(final_df['3Days_AccRF'] +  final_df['3D_RF_forecast'] >= final_df['Mean_Prec']), fbf_sug] = 'GO'
# final_df.loc[(final_df['3Days_AccRF'] +  final_df['5D_RF_forecast'] >= final_df['Mean_Prec']), fbf_sug] = 'READY'
# final_df.loc[(final_df['3Days_AccRF'] +  final_df['7D_RF_forecast'] >= final_df['Mean_Prec']), fbf_sug] = 'SET'
# final_df.loc[(final_df['3Days_AccRF'] +  final_df['10D_RF_forecast'] >= final_df['Mean_Prec']), fbf_sug] = 'MONITOR'
# final_df.loc[(final_df['3Days_AccRF'] +  final_df['10D_RF_forecast'] < final_df['Mean_Prec']), fbf_sug] = 'MONITOR'

#See notes for values
# final_df.loc[(final_df[rf_cat] == 'Light') & (final_df['NOV_%Normal_RF'] == 'Near Normal') &  (final_df['Flood_HAZ'] != 0), fbf_sug] = 'MONITOR'
# final_df.loc[(final_df[rf_cat] == 'Moderate') & (final_df['NOV_%Normal_RF'] == 'Near Normal') &  (final_df['Flood_HAZ'] != 0), fbf_sug] = 'MONITOR'
# final_df.loc[(final_df[rf_cat] == 'Intense') & (final_df['NOV_%Normal_RF'] == 'Near Normal') &  (final_df['Flood_HAZ'] != 0), fbf_sug] = 'READY'
# final_df.loc[(final_df[rf_cat] == 'Heavy') & (final_df['NOV_%Normal_RF'] == 'Near Normal') &  (final_df['Flood_HAZ'] != 0), fbf_sug] = 'SET'
# final_df.loc[(final_df[rf_cat] == 'Torrential') & (final_df['NOV_%Normal_RF'] == 'Above Normal') &  (final_df['Flood_HAZ'] != 0), fbf_sug] = 'GO'

#See notes for values
final_df.loc[(final_df["Flood_HAZ"] != 0) & (final_df['NOV_%Normal_RF'] == 'Below Normal') & (final_df[rf_cat] == 'Light') , fbf_sug] = 'MONITOR'
final_df.loc[(final_df["Flood_HAZ"] != 0) & (final_df['NOV_%Normal_RF'] == 'Near Normal') & (final_df[rf_cat] == 'Moderate'), fbf_sug] = 'MONITOR'
final_df.loc[(final_df["Flood_HAZ"] != 0) & (final_df['NOV_%Normal_RF'] == 'Near Normal') & (final_df[rf_cat] == 'Intense'), fbf_sug] = 'READY'
final_df.loc[(final_df["Flood_HAZ"] != 0) & (final_df['NOV_%Normal_RF'] == 'Above Normal') & (final_df[rf_cat] == 'Heavy'), fbf_sug] = 'SET'
final_df.loc[(final_df["Flood_HAZ"] != 0) & (final_df['NOV_%Normal_RF'] == 'Above Normal') & (final_df[rf_cat] == 'Torrential'), fbf_sug] = 'GO'

final_df = final_df.sort_values('3D_RF_forecast', ascending=[False])
final_df.dropna(subset=['3D_RF_forecast'], inplace=True)
final_df['Popn_at_Risk_FL'] = final_df['Popn_at_Risk_FL'].round()
cols = ['UID','Municipali','ADM4_EN','3Days_AccRF','10Days_AccRF','Mean_Prec','NOV_%Normal_RF','3D_RF_forecast','5D_RF_forecast',
        '7D_RF_forecast','10D_RF_forecast','Flood_HAZ','Popn_VUL_1','Popn_at_Risk_FL','No_of_Poor',rf_cat,fbf_sug]

cols_g = ['UID','Municipali','ADM4_EN','3Days_AccRF','10Days_AccRF','Mean_Prec','NOV_%Normal_RF','3D_RF_forecast','5D_RF_forecast',
        '7D_RF_forecast','10D_RF_forecast','Flood_HAZ','Popn_VUL_1','Popn_at_Risk_FL','No_of_Poor',rf_cat,fbf_sug, 'geometry']

df_clean_g = final_df[cols_g]
df_clean_g = df_clean_g.loc[df_clean_g[fbf_sug] != ""]

df_clean = final_df[cols]
df_clean = df_clean.loc[df_clean[fbf_sug] != ""]

df_clean.to_excel(alert_csv_path + alert_data_name)

b = gpd.GeoDataFrame(df_clean_g)
bbox = b.envelope
bbox_gdf = gpd.GeoDataFrame(gpd.GeoSeries(bbox), columns=['geometry'])

x = bbox.to_json()
m = folium.Map(tiles='cartodbpositron') #cartodbpositron #stamentoner #cartodbdark_matter
folium.GeoJson(x).add_to(m)
m.fit_bounds(m.get_bounds())
view = m.save(map_view_path + map_view_name)

job()
