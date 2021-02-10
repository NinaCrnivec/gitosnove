import epygram
import matplotlib.pyplot as plt
import numpy as np
import sys

epygram.init_env()

#---GLOBAL:

#Output variables:
var_P = "Total precipitation"
var_T = "CLSTEMPERATURE"
var_RH = "CLSHUMI.RELATIVE"
var_Td = "CLSDEWPOINTTEMPERATURE"
var_U = "CLSVENT.ZONAL"
var_V = "CLSVENT.MERIDIEN"
var_WS = "WINDSPEED"
var_WG = "10 metre wind gust"
var_C = "SURFNEBUL.TOTALE"
var_MSLP = "MSLPRESSURE"
var_Ed = "SURFRAYT SOLA DE"
var_En = "SURFFLU.RAY.SOLA"
var_G = "SURFGEOPOTENTIEL"

#Output variables on standard pressure levels:
varp_U = "WIND.U.PHY"
varp_V = "WIND.V.PHY"
varp_T = "TEMPERATUR"
varp_RH = "HUMI_RELAT"
varp_G = "GEOPOTENTI"
varp_W = "VITESSE_VE"

#Output file path:
fp_out="/scratch/ms/si/sinc/output_grib/"


def create_varset_lelamC(fc):
   var_arr = [var_T, var_Td, var_WS, var_G]
   if(fc>0):
      var_arr.remove(var_G)
      var_arr.insert(0, var_P)
      var_arr.extend([var_Ed, var_En])
   return var_arr

def create_varset_latlon(fc):
   var_arr = [var_T, var_RH, var_Td, var_U, var_V, var_WG, var_C, var_MSLP, var_G]
   if(fc>0):
      var_arr.remove(var_G)
      var_arr.insert(0, var_P)
   return var_arr

def create_generic_dict_P(fc):
   gen_dict_P = {
      'typeOfFirstFixedSurface': 1,
      'discipline': 0,
      'parameterCategory': 1,
      'parameterNumber': 52,
      'productDefinitionTemplateNumber': 8,
      'generatingProcessIdentifier' : 96,
      'hoursAfterDataCutoff' : 0,
      'minutesAfterDataCutoff' : 0,
      'typeOfStatisticalProcessing': 1,
      'level': 0,
      'indicatorOfUnitForTimeRange': 1,
      'typeOfTimeIncrement': 2,
      'forecastTime' : 0,
      'lengthOfTimeRange': fc,
      'stepUnits': 1,
      'startStep': 0,
      'endStep': fc,
   }
   return gen_dict_P

def create_generic_dict_E(paramNum, fc):
   gen_dict_E = {
      'discipline': 0,
      'parameterCategory': 4,
      'parameterNumber': paramNum,
      'typeOfFirstFixedSurface': 1,
      'typeOfStatisticalProcessing':1,
      'productDefinitionTemplateNumber': 8,
      'lengthOfTimeRange': fc,
      'startStep': 0,
      'endStep': fc,
      'stepUnits': 1,
   }
   return gen_dict_E

def create_generic_dict_UV_plev(paramNum, plev):
   gen_dict_UV_plev = {
      'discipline': 0,
      'parameterCategory': 2,
      'parameterNumber': paramNum,
      'productDefinitionTemplateNumber': 0,
      'level': plev,
      'indicatorOfTypeOfLevel':100,
      'typeOfFirstFixedSurface': 100
   }
   return gen_dict_UV_plev

class get_generic_dict:

   gen_dict_Td = {
      'typeOfFirstFixedSurface': 103,
      'level': 2,
      'discipline': 0,
      'parameterCategory': 0,
      'parameterNumber': 6,
      'productDefinitionTemplateNumber': 0,
   }

   gen_dict_WS = {
      'typeOfFirstFixedSurface': 103,
      'level': 10,
      'discipline': 0,
      'parameterCategory': 2,
      'parameterNumber': 1,
      'productDefinitionTemplateNumber': 0,
   }

   gen_dict_WG = {
      'typeOfFirstFixedSurface': 103,
      'discipline': 0,
      'parameterCategory': 2,
      'parameterNumber': 22,
      'productDefinitionTemplateNumber': 8,
      'typeOfStatisticalProcessing': 2,
      'lengthOfTimeRange': 1
   }

def reset_fid(arb_field, var_name, generic_dict):
   arb_field_fid = arb_field.fid
   arb_field_fid["FA"] = var_name
   arb_field_fid["generic"] = generic_dict
   arb_field.fid = arb_field_fid
   return arb_field

def set_var_fid(var, arb_field, fc):
   if var == var_P:
      generic_dict = create_generic_dict_P(fc)
   elif var == var_Td:
      generic_dict = get_generic_dict.gen_dict_Td
   elif var == var_WS:
      generic_dict = get_generic_dict.gen_dict_WS
   elif var == var_WG:
      generic_dict = get_generic_dict.gen_dict_WG
   elif var == var_Ed:
      generic_dict = create_generic_dict_E(7, fc)
   elif var == var_En:
      generic_dict = create_generic_dict_E(9, fc)
   else:
      pass

   if var==var_P or var==var_Td or var==var_WS or var==var_WG or var==var_Ed or var==var_En:
      arb_field = reset_fid(arb_field, var, generic_dict)

   return arb_field

def get_P(res):
   p1 = res.readfield("SURFPREC.EAU.CON")
   p2 = res.readfield("SURFPREC.NEI.CON")
   p3 = res.readfield("SURFPREC.EAU.GEC")
   p4 = res.readfield("SURFPREC.NEI.GEC")
   p = p1.deepcopy()
   p_data = p1.data + p2.data + p3.data + p4.data
   p.setdata(p_data)
   return p

def calc_Td(inp_T, inp_RH): # inp_T[K], inp_RH[/] between 0 and 1;
   #Define constants:
   T0  = 273.  # [K]
   es0 = 610.  # [Pa] 6.1 mbar=6.1 hPa = 610 Pa
   Rv  = 461.  # [J/kg/K]
   hi  = 2.5e6 # [J/kg] 2.5 MJ/kg
   #Calculate saturation water vapour pressure - es(T):
   es = es0*np.exp((hi/Rv)*(1/T0-1/inp_T))
   #Calculate partial water vapour pressure - e:
   e = inp_RH*es
   #Calculate Td:
   Td = (1/T0)-(Rv/hi)*np.log(e/es0)
   Td = 1.0/Td
   return Td

def get_dewpoint(T, RH):
   Td = T.deepcopy()
   Td_data = np.where(T.data > 0.0, calc_Td(T.data, RH.data), T.data)
   #Td_data = calc_Td(T.data, RH.data)
   Td.setdata(Td_data)
   return Td

def get_Td(res):
   T = res.readfield("CLSTEMPERATURE")
   RH = res.readfield("CLSHUMI.RELATIVE") #CLSHUMI.RELATIVE values between 0 and 1;
   Td = get_dewpoint(T, RH)
   return Td

def calc_wind_speed(u, v):
   wind_speed = u.deepcopy()
   wind_speed_data = np.sqrt((u.data)**2 + (v.data)**2)
   wind_speed.setdata(wind_speed_data)
   return wind_speed

def get_WS(res):
   u = res.readfield("CLSVENT.ZONAL")
   v = res.readfield("CLSVENT.MERIDIEN")
   windspeed = calc_wind_speed(u, v)
   return windspeed

def get_WG(res):
   u = res.readfield("CLSU.RAF.MOD.XFU")
   v = res.readfield("CLSV.RAF.MOD.XFU")
   windspeed = calc_wind_speed(u, v)
   return windspeed

def convert_to_percent(f):
   f_data = f.data*100.
   f.setdata(f_data)
   return f

def open_res(output_file):
   output = epygram.formats.resource(output_file, "w", fmt="GRIB")
   return output

def writefield_LJ_latlon(out_res, arb_field):
   out_res.writefield(arb_field, dict(centre=219))

def writefield_LJ_lelamC(out_res, arb_field):
   g_field = arb_field.geometry
   Ni = g_field.dimensions['X']
   Nj = g_field.dimensions['Y']
   Np = Ni*Nj
   out_res.writefield(arb_field, dict(centre=219, numberOfDataPoints=Np))

def read_field(var, res):
   if var == var_P:
      f = get_P(res)
   elif var == var_Td:
      f = get_Td(res)
   elif var == var_WS:
      f = get_WS(res)
   elif var == var_WG:
      f = get_WG(res)
   else:
      f = res.readfield(var)
   return f

def get_coordinates_lelamC(res):
   g = res.geometry
   numX_CI = g.dimensions["X_CIzone"]
   numY_CI = g.dimensions["Y_CIzone"]
   numX_Iwidth = g.dimensions["X_Iwidth"]
   numY_Iwidth = g.dimensions["Y_Iwidth"]
   ix_start = int(numX_Iwidth)
   ix_end   = int(numX_CI) - int(numX_Iwidth)
   iy_start = int(numY_Iwidth)
   iy_end   = int(numY_CI) - int(numY_Iwidth)
   coordinates = (ix_start, ix_end, iy_start, iy_end)
   return coordinates

class resampling_settings_latlon:
   borders = {
     "lonmin":  9.4,
     "latmin": 28.9,
     "lonmax": 45.0,
     "latmax": 52.5,
   }
   resolution_lon = 0.0300
   resolution_lat = 0.0225
   resampling_parameters = (borders, resolution_lon, resolution_lat)

def process_var_lelamC(r, fc):
   o = open_res(fp_out + "lelamC.grb")
   #Get parameters for the extraction of C zone out of CI zone on LELAM grid:
   (ix_start, ix_end, iy_start, iy_end) = get_coordinates_lelamC(r)
   var_arr_lelamC = create_varset_lelamC(fc)
   for var in var_arr_lelamC:
      print(var)
      f = read_field(var, r)
      f.select_subzone(subzone="CI")
      f_lelamC = f.extract_subarray(ix_start, ix_end, iy_start, iy_end)
      f_lelamC = set_var_fid(var, f_lelamC, fc)
      if var == var_G:
         og = open_res(fp_out + "alsmsw_lela_geop.grb2")
         writefield_LJ_lelamC(og, f_lelamC)
      else:
         writefield_LJ_lelamC(o, f_lelamC)
   return


def process_var_latlon(r, fc):
   o  = open_res(fp_out + "latlon.grb")
   #Set parameters for resampling on regular lat-lon grid:
   (borders, resolution_lon, resolution_lat) = resampling_settings_latlon.resampling_parameters
   var_arr_latlon = create_varset_latlon(fc)
   for var in var_arr_latlon:
      print(var)
      if var == var_Td:
         f = get_dewpoint(T_latlon, RH_latlon) #Calculate Td from resampled T and RH [RH between 0 and 1]
      else:
         f = read_field(var, r)
      f_latlon = f.resample_on_regularll_mod(borders, resolution_lon, resolution_lat)
      if var == var_T:  T_latlon = f_latlon
      if var == var_RH: RH_latlon = f_latlon
      if var == var_C or var == var_RH:
         f_latlon = convert_to_percent(f_latlon)
      f_latlon = set_var_fid(var, f_latlon, fc)
      if var == var_G:
         og = open_res(fp_out + "alsmsw_lalo_geop.grb2")
         writefield_LJ_latlon(og, f_latlon)
      else:
         writefield_LJ_latlon(o, f_latlon)


def process_var_latlon_plev(r, fc):
   #Set parameters for resampling on regular lat-lon grid:
   (borders, resolution_lon, resolution_lat) = resampling_settings_latlon.resampling_parameters
   plev_arr = ["20000","25000","30000","40000","50000","60000","70000","80000","85000","90000","92500","95000","00000"]
   #plev_arr = ["50000"] #test version
   varp_arr = [varp_U, varp_V, varp_T, varp_RH, varp_G, varp_W]
   for plev in plev_arr:
      print("plev = ", plev)
      op_name = fp_out + "fp_latlon_P" + plev
      op = open_res(op_name + ".grb")
      for var in varp_arr:
         print("var = ", var)
         if var == varp_U:
            opu = open_res(op_name + "_u.grb")
         elif var == varp_V:
            opv = open_res(op_name + "_v.grb")
         else:
             pass
         varname = "P" + plev + var
         f = r.readfield(varname)
         if var == varp_RH:
            f = convert_to_percent(f)
         f_latlon = f.resample_on_regularll_mod(borders, resolution_lon, resolution_lat)
         if var == varp_U:
            generic_dict = create_generic_dict_UV_plev(2, plev)
            f_latlon = reset_fid(f_latlon, var, generic_dict)
            writefield_LJ_latlon(opu, f_latlon)
         elif var == varp_V:
            generic_dict = create_generic_dict_UV_plev(3, plev)
            f_latlon = reset_fid(f_latlon, var, generic_dict)
            writefield_LJ_latlon(opv, f_latlon)
         else:
            writefield_LJ_latlon(op, f_latlon)


def check_input(input_file, script_file):
   if(input_file==script_file):
      print("Please enter input file name!")
      sys.exit()
   try:
      file = open(input_file, "r")
   except IOError:
      print("There was an error reading file!")
      sys.exit()


def main():

   #inp_file = "/scratch/ms/si/sism/rundir/as2020112200_seemhews_oper_prod/seemhews_oper/prod/forecast/integration/PFSMWSLELA+0072"

   #Input file:   
   inp_file = sys.argv[-1] #takes the last argument from the command line (returns the scriptfile if nothing is entered)  
   script_file = sys.argv[0]
   check_input(inp_file, script_file)

   #Open resource:
   r = epygram.formats.resource(inp_file, "r")
   fc = r.validity.term(fmt="IntHours")

   #Process variables on lelamC:
   process_var_lelamC(r, fc)

   #Process variables on latlon:
   process_var_latlon(r, fc)

   #Process variables on latlon on pressure levels:
   process_var_latlon_plev(r, fc)



if __name__ == "__main__":
    main()

















