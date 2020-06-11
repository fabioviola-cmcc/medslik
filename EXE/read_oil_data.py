import sys

path="oilbase.txt"
data=open(path).readlines()

#oil_list=""
#for i in range(0,len(data)):
# element=data[i].split('  ')
# name=element[0]
# oil_list=oil_list+name+"\n"

#doc0=file("oil_list.txt",'w')
#doc0.writelines(oil_list)
#doc0.close()

option=sys.argv[1]

if option=="NAME":
 oil_name=sys.argv[2]
#oil_name="Aboozar"
 for i in range(0,len(data)):
  element_=data[i].split('  ')
  name=element_[0]
  element=data[i].split(' ')
  if name==oil_name:
   n=len(element)-1
   dens_oil=element[n][0:5]
   api=element[n][5:10]
   temp=element[n][10:15]
   visc=element[n][18:23]
   res_dens_oil=element[n][23:28]
   res_perc_oil=element[n][28:33]
   vap_press=element[n][33:38]
   oil_file=name+"\n" 
   oil_file=oil_file+api+"          API of Oil\n"
   oil_file=oil_file+dens_oil+"          Density of Oil\n"
   oil_file=oil_file+res_dens_oil+"          Residual Density of Oil\n"
   oil_file=oil_file+res_perc_oil+"          Residual Percent of Oil\n"
   oil_file=oil_file+visc+"          Viscosity of Oil\n"
   oil_file=oil_file+temp+"          Temperature at which Viscosity determined\n"
   oil_file=oil_file+vap_press+"          Vapour Pressure of Oil (bar)\n"
   doc=file("oil_file.txt",'w')
   doc.writelines(oil_file)
   doc.close()

if option=="API":
 api_number=sys.argv[2]
#api_number="20"
 for i in range(0,223):
  element=data[i].split(' ')
  name=element[0]
  n=len(element)-1
  api=element[n][5:10]
  api_int=str(round(float(api)))[0:2]
  if api_int==str(round(float(api_number)))[0:2]:
#   print api_number
   dens_oil=element[n][0:5]
   api=element[n][5:10]
   temp=element[n][10:15]
   visc=element[n][18:23]
   res_dens_oil=element[n][23:28]
   res_perc_oil=element[n][28:33]
   vap_press=element[n][33:38]
   oil_file=name+"\n" 
   oil_file=oil_file+api+"          API of Oil\n"
   oil_file=oil_file+dens_oil+"          Density of Oil\n"
   oil_file=oil_file+res_dens_oil+"          Residual Density of Oil\n"
   oil_file=oil_file+res_perc_oil+"          Residual Percent of Oil\n"
   oil_file=oil_file+visc+"          Viscosity of Oil\n"
   oil_file=oil_file+temp+"          Temperature at which Viscosity determined\n"
   oil_file=oil_file+vap_press+"          Vapour Pressure of Oil (bar)\n"
   doc=file("oil_file.txt",'w')
   doc.writelines(oil_file)
   doc.close()
