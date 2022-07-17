import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import statistics
from scipy import stats
import math

import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import locale

DICT_PLOTMODE={"00":"T;P_iz","01":"T;U_iz","02":"T;I_iz","10":"T;η","11":"T;η","12":"T;η","20":"E;P_iz","21":"E;U_iz","22":"E;I_iz","30":"E;η","31":"E;η","32":"E;η","40":"E;T","41":"E;T","42":"E;T"}
DICT_TYPE_OF_PV_PANEL={'REF':'REFERENTNI','HLA':'HLAĐENI'}
DICT_DATE_MONTH={"1":"Siječanj","2":"Veljača","3":"Ožujak","4":"Travanj","5":"Svibanj","6":"Lipanj","7":"Srpanj","8":"Kolovoz","9":"Rujan","10":"Listopad","11":"Studeni","12":"Prosinac"}
DICT_PLOTLABEL={"00":"izlazna snaga","01":"izlazni napon","02":"izlazna struja","10":"učinkovitost","11":"učinkovitost","12":"učinkovitost"}

AXIS_VALUE_PRECISION=0

#USER DEFINED:

LABEL_FONT_SIZE=12

SHOW_LSR_LINE = 1
SHOW_ONLY_LSR_LINE=0 #dont use in combination with ADJUSTABLE_Y_AXIS==1

MULTIPLE_DAYS_GRAPH = 0#Show more graphs of each day in one plot
MULTIPLE_MONTHS_GRAPH = 0 #Show more graphs of each month in one plot

ADJUSTABLE_X_AXIS = 1 #for solving saturation problems of DATE type values
TIME_SCALE=4 #interval between dates in x axis; must be ADJUSTABLE_X_AXIS=1
ADJUSTABLE_Y_AXIS = 1 #0->Python-AUTO mode, 1-> y values from 0 to max

REMOVE_FAULT_DATA=1 #Removes measurings by large deviations

MONTHS=[8,8] #Range of months
DAYS=[21,21] #Range of days

TYPE_OF_PV_PANEL='REF' #REF->referentni FN modul, HLA->hladjeni FN modul

SET_IRRADIANCE = 800#in W/m2 ; use -1 for acquiring data for all values of irradiance PLOT_MODE(<5)
SET_IRRADIANCE_VALUES=[600,1100,10] #[lowerBound,upperBound,step] USE ONLY FOR CALCULATION OF TEMPERATURE COEFFICIENTS

IRRADIANCE_LOWER_BOUND_FACTOR = 0.995
IRRADIANCE_UPPER_BOUND_FACTOR = 1.005

SET_TEMPERATURE = -1 #in oC; use -1 for acquiring data for all values of temperature (for example temperature coefficients)

TEMPERATURE_LOWER_BOUND_FACTOR = 0.995
TEMPERATURE_UPPER_BOUND_FACTOR = 1.005

PLOT_MODE=1 #0->plot electrical parameters (temperature) (voltage,current,output_power), 1->plot efficiency(temperature)
            #2-> plot electrical parameters (irradiance) (voltage,current,output_power), 3->plot efficiency(irradiance), 4->temperature-irradiance (use for irradiance SET_IRRADIANCE_VALUES)
ADDITIONAL_MODE =0#0->None, 1-> 3D plot, 2->temperature efficiency plot with irradiance as additional parameter (2D plot),
                  #3->plot max.efficiency for date & independent and dependent histograms of temperature and irradiance, 4-> time plots (use TIME_MODE)
TIME_MODE=0#0->irradiance, 1->temperature, 2->electrical parameters, 3->efficiency
ELECTRICAL_PARAM_MODE=0 #0->output_power, 1->voltage, 2->current
TEMPERATURE_MODE=1 #0-> UR REF, 1-> C REF, 2-> DL REF, 3-> AVG REF

AREA_OF_PV_PANEL = 1.62688 ##m^2
SHOW_LOG=0#0-> Not show log data in shell, 1-> Show log data in shell
SHOW_STATISTICS = 1#Show statistical data of parameters, LSR equation, table view of best efficiencies of days...
SHOW_DATA_FOR_EXPORT = 1#exports print of the table with date-temperature-irradiance-max efficiency for data analysis


def get_data_with_set_irradiance(data_path, input_irradiance,factor_lower_bound,factor_higher_bound):
    _result=[]

    if(float(input_irradiance)==-1.0):
        factor_lower_bound=float(0)
        factor_higher_bound=float('-inf')
    
    try:
        with open(data_path) as _csv_file:
            _csv_reader=csv.reader(_csv_file,delimiter=',')
            _line_count=0
            for row in _csv_reader:
                if(_line_count>=1):
                    _date_comp=row[0].split("-") ##yyyy-mm-dd -> yyyy, mm, dd because of US style and zero before date <=9
                    _time=_date_comp[2]+"."+_date_comp[1]+"."+_date_comp[0]+"#"+(row[1])[1:len(row[1])] ##6.8.2021#05:35:00
                    if((_time.split("#")[1])[0:5]!='00:00'):
                        if(float(row[2])>=factor_lower_bound*input_irradiance)and(float(row[2])<=factor_higher_bound*input_irradiance):
                            _result.append(row[2]) #need to convert to float back
                            _result.append(_time)
                    else: break
                _line_count+=1
            #print("Number of data in Apogee.csv: "+str(len(result)))
    except Exception as e:
        if(SHOW_LOG==1): print(e)

    return _result

def get_electrical_data_avg(data_path, irradianceData, mode):
    _result=[]
    _modeText=""
    _line_count=0
    
    try:
        with open(data_path) as _csv_file:
            _csv_reader = csv.reader(_csv_file, delimiter=';') ##ATTENTION for file M0101\REF.csv delimiter is ';', for file M0815\REF.csv and other files delimiter is '\t'
            _line_count=0
            for row in _csv_reader:
                if(_line_count>0):
                    _time_comp=row[1].split("-")
                    _time=(_time_comp[2])[0:2]+"."+_time_comp[1]+"."+_time_comp[0]+"#"+(_time_comp[2])[3:11]
                    if((_time.split("#")[1])[0:5]!='00:00'):
                        for i in range(1,len(irradianceData),2):
                            if(str(_time)[0:16]==str(irradianceData[i])[0:16]):
                                if(mode==0):#power
                                    _result.append(str(float(row[3])))
                                    _result.append(_time)
                                elif(mode==1):  #voltage                      
                                    _result.append(str(-float(row[6])))
                                    _result.append(_time)
                                elif(mode==2): #current
                                    _result.append(str(-float(row[12])/1000.0))
##                                      _result.append(str(float(row[3])/-float(row[6])))
                                    _result.append(_time)
                                else:
                                    print("Mod nije dostupan!")
                                    _result.append(_time)
                    else: break
                       
                _line_count+=1
                
    except Exception as e:
        if(SHOW_LOG==1): print(e)

    if(mode==0):
        _modeText="Prosječna proizvedena snaga FN modula,P [W]"
    elif(mode==1):
        _modeText="Prosječan iznos napona FN modula,U [V]"
    elif(mode==2):
        _modeText="Prosječna jakost struje FN modula,I [A]"
    else:
        print("Mod nije dostupan!")

    return _result,_modeText

def get_input_power(input_irradiance_data):
    _result=[]
    _textMode=""

    for i in range(0,len(input_irradiance_data),2):
        _result.append(str(float(input_irradiance_data[i])*float(AREA_OF_PV_PANEL)))
        _result.append(input_irradiance_data[i+1])

    _textMode="Ulazna snaga FN modula [W]"

    return _result,_textMode

def get_efficiency_irradiance(month,day,input_irradiance):
    _result=[]
    _output_power=[]
    _input_power=[]
    _modeText=""

    _output_power=get_electrical_data_avg(get_data_path(month,day,0),input_irradiance,0)[0]    
    _input_power=get_input_power(input_irradiance)[0]

    try:
        for i in range(1,len(_input_power),2):
            for j in range(1,len(_output_power),2):
                if(str(_input_power[i])[0:16]==str(_output_power[j])[0:16]):
                    _eff_val=0.0
                    if(float(_input_power[i-1])!=0.0):
                        _eff_val=float(_output_power[j-1])/float(_input_power[i-1])
                        _result.append(str(_eff_val))
                        _result.append(_input_power[i])

    except Exception as e:
        if(SHOW_LOG==1): print(e)

    _modeText="Učinkovitost FN modula,η"

    return _result,_modeText


def get_data_with_set_temperature(data_path, mode, input_temperature,factor_lower_bound,factor_higher_bound):
    _result=[]

    _REF_FRONT_UR=[]
    _REF_FRONT_C=[]
    _REF_FRONT_DL=[]

    _HLA_FRONT_UR=[]
    _HLA_FRONT_C=[]
    _HLA_FRONT_DL=[]

    if(float(input_temperature)==-1.0):
        factor_lower_bound=float(0)
        factor_higher_bound=float('-inf')
    
    try:
        _currentFile=open(data_path,'r')
        _lines=_currentFile.readlines()
        _counter=0
        for line in _lines:
            if(_counter>1):
                _components_of_line=line.split("\t")
                _date=_components_of_line[0].replace("/",".")
                _time_comp=_components_of_line[1]
                _time=_date+"#"+_time_comp
                if((_time.split("#")[1])[0:5]!='00:00'):
                #print(time)
                    if(mode==0):
                        if(TYPE_OF_PV_PANEL=='REF' and (float(_components_of_line[5])>=factor_lower_bound*input_temperature)and(float(_components_of_line[5])<=factor_higher_bound*input_temperature)):
                            _result.append(_components_of_line[5])
                            _result.append(_time)
                        elif(TYPE_OF_PV_PANEL=='HLA' and (float(_components_of_line[2])>=factor_lower_bound*input_temperature)and(float(_components_of_line[2])<=factor_higher_bound*input_temperature)):
                            _result.append(_components_of_line[2])
                            _result.append(_time)
                    elif(mode==1):
                        if(TYPE_OF_PV_PANEL=='REF'and (float(_components_of_line[6])>=factor_lower_bound*input_temperature)and(float(_components_of_line[6])<=factor_higher_bound*input_temperature)):
                            _result.append(_components_of_line[6])
                            _result.append(_time)
                        elif(TYPE_OF_PV_PANEL=='HLA'and (float(_components_of_line[3])>=factor_lower_bound*input_temperature)and(float(_components_of_line[3])<=factor_higher_bound*input_temperature)):
                            _result.append(_components_of_line[3])
                            _result.append(_time)
                    elif(mode==2):
                        if(TYPE_OF_PV_PANEL=='REF'and (float(_components_of_line[7])>=factor_lower_bound*input_temperature)and(float(_components_of_line[7])<=factor_higher_bound*input_temperature)):
                            _result.append(_components_of_line[7])
                            _result.append(_time)
                        elif(TYPE_OF_PV_PANEL=='HLA'and (float(_components_of_line[4])>=factor_lower_bound*input_temperature)and(float(_components_of_line[4])<=factor_higher_bound*input_temperature)):
                            _result.append(_components_of_line[4])
                            _result.append(_time)                    
                    elif(mode==3):
                        if(TYPE_OF_PV_PANEL=='REF'and (float(_components_of_line[5])>=factor_lower_bound*input_temperature)and(float(_components_of_line[5])<=factor_higher_bound*input_temperature)and (float(_components_of_line[6])>=factor_lower_bound*input_temperature)and(float(_components_of_line[6])<=factor_higher_bound*input_temperature)and (float(_components_of_line[7])>=factor_lower_bound*input_temperature)and(float(_components_of_line[7])<=factor_higher_bound*input_temperature)):
                            _REF_FRONT_UR.append(float(_components_of_line[5]))
                            _REF_FRONT_C.append(float(_components_of_line[6]))
                            _REF_FRONT_DL.append(float(_components_of_line[7]))
                            _result.append(str((_REF_FRONT_UR[_counter-2]+_REF_FRONT_C[_counter-2]+_REF_FRONT_DL[_counter-2])/3))
                            _result.append(_time)
                        elif(TYPE_OF_PV_PANEL=='HLA'and (float(_components_of_line[5])>=factor_lower_bound*input_temperature)and(float(_components_of_line[5])<=factor_higher_bound*input_temperature)and (float(_components_of_line[6])>=factor_lower_bound*input_temperature)and(float(_components_of_line[6])<=factor_higher_bound*input_temperature)and (float(_components_of_line[7])>=factor_lower_bound*input_temperature)and(float(_components_of_line[7])<=factor_higher_bound*input_temperature)):
                            _HLA_FRONT_UR.append(float(_components_of_line[2]))
                            _HLA_FRONT_C.append(float(_components_of_line[3]))
                            _HLA_FRONT_DL.append(float(_components_of_line[4]))
                            _result.append(str((_HLA_FRONT_UR[_counter-2]+_HLA_FRONT_C[_counter-2]+_HLA_FRONT_DL[_counter-2])/3))
                            _result.append(_time)
                    else:
                        _result.append(0)
                        _result.append(_time)
                else: break
                
            _counter+=1
            
    except Exception as e:
        if(SHOW_LOG==1): print(e)
        
    if(mode==0):
        if(TYPE_OF_PV_PANEL=='REF'):
            _modeText="REF_FRONT_UR=Temperatura na gornje-desnom kutu površine FN modula,T [°C]"
        else: _modeText="HLA_FRONT_UR=Temperatura na gornje-desnom kutu površine FN modula,T [°C]"
    elif(mode==1):
        if(TYPE_OF_PV_PANEL=='REF'):
            _modeText="REF_FRONT_C=Temperatura u sredini površine FN modula,T [°C]"
        else: _modeText="HLA_FRONT_C=Temperatura u sredini površine FN modula,T [°C]"
    elif(mode==2):
        if(TYPE_OF_PV_PANEL=='REF'):
            _modeText="REF_FRONT_DL=Temperatura na donje-lijevom kutu površine FN modula,T [°C]"
        else: _modeText="HLA_FRONT_DL=Temperatura na donje-lijevom kutu površine FN modula,T [°C]"
    elif(mode==3):
        if(TYPE_OF_PV_PANEL=='REF'):
            _modeText="REF_FRONT_AVG=Prosječna temperatura površine FN modula,T [°C]"
        else: _modeText="HLA_FRONT_AVG=Prosječna temperatura površine FN modula,T [°C]"
    else: _modeText="Nije dostupno"
            
    return _result, _modeText

def get_data_path(month,day,mode):
    _data_path=""
    
    if (month>=1 and month<=9)and(day>=1 and day<=9):
        _data_path='M0'+str(month)+'0'+str(day)
    elif day>=1 and day<=9:
        _data_path='M'+str(month)+'0'+str(day)
    elif month>=1 and month<=9:
        _data_path='M0'+str(month)+str(day)
    else:
        _data_path='M'+str(month)+str(day)

    if(mode==0):
        _data_path=_data_path+'\\'+TYPE_OF_PV_PANEL+ '.csv'
    elif(mode==1):
        _data_path=_data_path+'\paneli.txt' ##ATTENTION some folders Mmmdd dont have this file
    elif(mode==2):
        _data_path=_data_path+'\Apogee.csv'
        
    return _data_path

def match_XY_plot_data(arrayX,arrayY):
    _resultX_value=[]
    _resultY_value=[]
    
    try:
        if ((arrayX) and (arrayY)):
            for i in range(1,len(arrayX),2):
                for j in range(1,len(arrayY),2):
                    if(str(arrayX[i])[0:16]==str(arrayY[j])[0:16]):
                        if(SHOW_LOG==1):
                            print("X_TIME: "+str(arrayX[i])[0:16]+" X_VALUE: "+str(arrayX[i-1])+" Y_TIME: "+str(arrayY[j])[0:16]+" Y_VALUE: "+str(arrayY[j-1]))
                        _resultX_value.append(float(arrayX[i-1]))
                        _resultY_value.append(float(arrayY[j-1]))
        else:
            if(SHOW_LOG==1): print("No data found")
            return _resultX_value,_resultY_value

    except Exception as e:
        if(SHOW_LOG==1): print(e)

    if(not(len(_resultX_value) == len(_resultY_value))):
        _resultX_value=[]
        _resultY_value=[]

    return _resultX_value,_resultY_value

def match_XYZ_plot_data(arrayX,arrayY,arrayZ):
    _resultX_value=[]
    _resultY_value=[]
    _resultZ_value=[]

    _resultY_data=[]

    try:

        if((arrayX) and (arrayY) and (arrayZ)):
            for i in range(1,len(arrayX),2):
                for j in range(1,len(arrayY),2):
                    if(str(arrayX[i])[0:16]==str(arrayY[j])[0:16]):
                        if(SHOW_LOG==1):
                            print("X_TIME: "+str(arrayX[i])[0:16]+" X_VALUE: "+str(arrayX[i-1])+" Y_TIME: "+str(arrayY[j])[0:16]+" Y_VALUE: "+str(arrayY[j-1]))
                        if(not(float(arrayX[i-1]) in _resultX_value) and not(arrayY[j-1] in _resultY_data)):
                            _resultX_value.append(float(arrayX[i-1]))                        
                            _resultY_data.append(arrayY[j-1])
                            _resultY_data.append(arrayY[j])

            for j in range(1,len(_resultY_data),2):
                for k in range(1,len(arrayZ),2):
                    if(str(_resultY_data[j])[0:16]==str(arrayZ[k])[0:16]):
                        if(SHOW_LOG==1):
                            print("Y_TIME: "+str(_resultY_data[j])[0:16]+" Y_VALUE: "+str(_resultY_data[j-1])+" Z_TIME: "+str(arrayZ[k])[0:16]+" Z_VALUE: "+str(arrayZ[k-1]))
                        if(not(float(_resultY_data[j-1]) in _resultY_value) and not(float(arrayZ[k-1]) in _resultZ_value)):
                            _resultY_value.append(float(_resultY_data[j-1]))
                            _resultZ_value.append(float(arrayZ[k-1]))

        else:
            if(SHOW_LOG==1): print("No data found")
            return _resultX_value,_resultY_value, _resultZ_value
                                             
    except Exception as e:
        if(SHOW_LOG==1): print(e)

    if(not(len(_resultX_value) == len(_resultY_value) and len(_resultY_value) == len(_resultZ_value))):
        _resultX_value=[]
        _resultY_value=[]
        _resultZ_value=[]
    else:
        if(SHOW_LOG==1):
            for i in range(len(_resultX_value)):
                print("X_VALUE: "+str(_resultX_value[i])+" Y_VALUE: "+str(_resultY_value[i])+" Z_VALUE: "+str(_resultZ_value[i]))

    return _resultX_value,_resultY_value,_resultZ_value

def date_exists(month,day):
    if(month==2):
        if(day>=1 and day<=28):
            return True
        else: return False
    elif(month<=7):
        if(month%2==1 and day>=1 and day<=31):
            return True
        elif(month%2==0 and day>=1 and day<=30):
            return True
        else: return False
    elif(month>=8):
        if(month%2==0 and day>=1 and day<=31):
            return True
        elif(month%2==1 and day>=1 and day<=30):
            return True
        else: return False
    else: return False

def fit_LSR(x,y,plotLine,plotOnlyLine,multipleMode,legendText):

    _coeffs=[0.0,0.0,0.0,0.0,0.0]
    _x=np.array(x)

    if((x) and (y)):
        _coeffs=stats.linregress(x,y)

        if(not(multipleMode)):
            if(plotOnlyLine==True):
                plt.plot(_x,_coeffs[0]*_x + _coeffs[1],color='r',label=legendText)
            else:    
                if(plotLine==True):
                    plt.plot(_x,_coeffs[0]*_x + _coeffs[1],color='r')        
        else:
            if(plotOnlyLine==True):
                plt.plot(_x,_coeffs[0]*_x + _coeffs[1],label=legendText)
            else:    
                if(plotLine==True):
                    plt.plot(_x,_coeffs[0]*_x + _coeffs[1])

    return _coeffs

def get_statistics(x,y):
    _statistics_data=[] #0,3->mean(x,y), 1,4->variance(x,y), 2,5->standard deviation(x,y)
    _statistics_data.append(statistics.mean(x))
    _statistics_data.append(statistics.variance(x))
    _statistics_data.append(statistics.stdev(x))
    _statistics_data.append(statistics.mean(y))
    _statistics_data.append(statistics.variance(y))
    _statistics_data.append(statistics.stdev(y))

    print('Prosjek od '+DICT_PLOTMODE[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].split(";")[0]+": "+str(_statistics_data[0])+'\n'+'Varijanca od '+DICT_PLOTMODE[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].split(";")[0]+': '+str(_statistics_data[1])+'\n'+'Standardna devijacija od '+DICT_PLOTMODE[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].split(";")[0]+': '+str(_statistics_data[2]))
    print('Prosjek od '+DICT_PLOTMODE[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].split(";")[1]+": "+str(_statistics_data[3])+'\n'+'Varijanca od '+DICT_PLOTMODE[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].split(";")[1]+': '+str(_statistics_data[4])+'\n'+'Standardna devijacija od '+DICT_PLOTMODE[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].split(";")[1]+': '+str(_statistics_data[5]))

    if(PLOT_MODE==0 or PLOT_MODE==1):
        print('Ozračenje E: '+str(SET_IRRADIANCE)+"\nPreciznost ozračenja E: "+str(IRRADIANCE_LOWER_BOUND_FACTOR)+"-"+str(IRRADIANCE_UPPER_BOUND_FACTOR))
    else:
        print('Interval ozračenja E: '+str(SET_IRRADIANCE_VALUES[0])+"-"+str(SET_IRRADIANCE_VALUES[1])+" s korakom: "+str(SET_IRRADIANCE_VALUES[2])+"\nPreciznost ozračenja E: "+str(IRRADIANCE_LOWER_BOUND_FACTOR)+"-"+str(IRRADIANCE_UPPER_BOUND_FACTOR))

    return _statistics_data
    
def remove_fault_data(x,y):
    try:
        _q1,_q2,_q3=statistics.quantiles(y,n=4)
        _irq=_q3-_q1
        _min_val=_q1-(1.5*_irq)
        _max_val=_q3+(1.5*_irq)

        if(SHOW_LOG==1):
            print("First quartile Q1: "+str((_min_val))+" Third quartile  Q3: "+str((_max_val)))
            
        _fault_indexes=[]
      
        for i in range(len(y)):
            if(float(y[i])<float(_min_val) or float(y[i])>float(_max_val) or math.isnan(float(y[i]))):
                _fault_indexes.append(i)
        
        for index in sorted(_fault_indexes, reverse=True):
            del x[index]
            del y[index]

    except Exception as e:
        if(SHOW_LOG==1): print(e)


def remove_fault_data_XYZ(x,y,z):
    try:
        _q1,_q2,_q3=statistics.quantiles(z,n=4)
        _irq=_q3-_q1
        _min_val=_q1-(1.5*_irq)
        _max_val=_q3+(1.5*_irq)

        if(SHOW_LOG==1):
            print("First quartile Q1: "+str((_min_val))+" Third quartile  Q3: "+str((_max_val)))

        _fault_indexes=[]
      
        for i in range(len(z)):
            if(float(z[i])<float(_min_val) or float(z[i])>float(_max_val) or math.isnan(float(z[i]))):
                _fault_indexes.append(i)
        
        for index in sorted(_fault_indexes, reverse=True):
            del x[index]
            del y[index]
            del z[index]

    except Exception as e:
        if(SHOW_LOG==1): print(e)

def show_fault_data_XYZDate(x,y,z,date):
    try:
        _q1,_q2,_q3=statistics.quantiles(z,n=4)
        _irq=_q3-_q1
        _min_val=_q1-(1.5*_irq)
        _max_val=_q3+(1.5*_irq)

        if(SHOW_LOG==1):
            print("First quartile Q1: "+str((_min_val))+" Third quartile  Q3: "+str((_max_val)))

        _fault_indexes=[]
      
        for i in range(len(z)):
            if(float(z[i])<float(_min_val) or float(z[i])>float(_max_val) or math.isnan(float(z[i]))):
                _fault_indexes.append(i)
        
        for index in sorted(_fault_indexes, reverse=False):
            print("Datum:"+date[index]+"  Temperatura [oC]:"+str(x[index])+"  Ozračenje [W/m2]:"+str(y[index])+"  Učinkovitost:"+str(z[index]))

    except Exception as e:
        if(SHOW_LOG==1): print(e)


def plot_mode_0(month, day, irradiance_val, irradiance_lower_bound_factor, irradiance_upper_bound_factor, temperature_mode, temperature_val, temperature_lower_bound_factor, temperature_upper_bound_factor, electrical_param_mode):

    _xx=[]
    _yy=[]
    _irradiance=[]
    _partialModeText=""
    _details=""
    _xLabel=""
    _yLabel=""
    
    _irradiance=get_data_with_set_irradiance(get_data_path(month, day,2), irradiance_val,irradiance_lower_bound_factor, irradiance_upper_bound_factor)
    _xx,_partialModeText=get_data_with_set_temperature(get_data_path(month, day,1), temperature_mode, temperature_val, temperature_lower_bound_factor, temperature_upper_bound_factor)
    _xLabel=_partialModeText.split("=")[1]
    _yy,_partialModeText=get_electrical_data_avg(get_data_path(month, day,0), _irradiance, electrical_param_mode)
    _yLabel=_partialModeText
    _details=("Ozračenje: "+str(irradiance_val)+" [W/m2]")+"     "+("Preciznost ozračenja: ["+str(IRRADIANCE_LOWER_BOUND_FACTOR)+","+str(IRRADIANCE_UPPER_BOUND_FACTOR)+"]")

    global AXIS_VALUE_PRECISION
    AXIS_VALUE_PRECISION=2

    return _xx,_yy,_xLabel,_yLabel,_details

def plot_mode_1(month, day, irradiance_val, irradiance_lower_bound_factor, irradiance_upper_bound_factor, temperature_mode, temperature_val, temperature_lower_bound_factor, temperature_upper_bound_factor):

    _xx=[]
    _yy=[]
    _irradiance=[]
    _partialModeText=""
    _details=""
    _xLabel=""
    _yLabel=""

    _irradiance=get_data_with_set_irradiance(get_data_path(month, day,2), irradiance_val, irradiance_lower_bound_factor, irradiance_upper_bound_factor)
    _xx,_partialModeText=get_data_with_set_temperature(get_data_path(month, day,1),temperature_mode, temperature_val, temperature_lower_bound_factor, temperature_upper_bound_factor)
    _xLabel=_partialModeText.split("=")[1]
    _yy,_partialModeText=get_efficiency_irradiance(month, day,_irradiance)
    _yLabel= "Učinkovitost FN modula,η"
    _details=("Ozračenje: "+str(irradiance_val)+" [W/m2]")+"     "+("Preciznost ozračenja: ["+str(IRRADIANCE_LOWER_BOUND_FACTOR)+","+str(IRRADIANCE_UPPER_BOUND_FACTOR)+"]")

    global AXIS_VALUE_PRECISION
    AXIS_VALUE_PRECISION=4

    return _xx,_yy,_xLabel,_yLabel,_details


def plot_mode_2(month, day, irradiance_lower_bound_factor, irradiance_upper_bound_factor, irradiance_values_start, irradiance_values_end, irradiance_values_step,electrical_param_mode):

    _xx=[]
    _yy=[]
    _x=[]
    _y=[]
    _partialModeText=""
    _modeTextString=""
    _details=""
    _xLabel=""
    _yLabel=""
    
    for irradiance_value in range(int(irradiance_values_start),int(irradiance_values_end),int(irradiance_values_step)):
        _x=get_data_with_set_irradiance(get_data_path(month, day,2), float(irradiance_value), irradiance_lower_bound_factor, irradiance_upper_bound_factor)
        _xLabel="Ozračenje FN modula,E [W/m2]"
        _y,_yLabel=get_electrical_data_avg(get_data_path(month, day,0), _x, electrical_param_mode)
        _details=("Ozračenje: ["+str(irradiance_values_start)+"-"+str(irradiance_values_end)+"] [W/m2]")+"     "+("Preciznost ozračenja: ["+str(IRRADIANCE_LOWER_BOUND_FACTOR)+","+str(IRRADIANCE_UPPER_BOUND_FACTOR)+"]")

        for j in _x:
            _xx.append(j)
        for j in _y:
            _yy.append(j)
            
    global AXIS_VALUE_PRECISION
    AXIS_VALUE_PRECISION=2
    
    return _xx,_yy,_xLabel,_yLabel,_details

def plot_mode_3(month, day, irradiance_lower_bound_factor, irradiance_upper_bound_factor, irradiance_values_start, irradiance_values_end, irradiance_values_step):

    _xx=[]
    _yy=[]
    _x=[]
    _y=[]
    _partialModeText=""
    _modeTextString=""
    _details=""
    _xLabel=""
    _yLabel=""
    
    for irradiance_value in range(int(irradiance_values_start),int(irradiance_values_end),int(irradiance_values_step)):
        _x=get_data_with_set_irradiance(get_data_path(month, day,2), float(irradiance_value),irradiance_lower_bound_factor, irradiance_upper_bound_factor)                    
        _y,_partialModeText=get_efficiency_irradiance(month, day,_x)
        _details=("Ozračenje: ["+str(irradiance_values_start)+"-"+str(irradiance_values_end)+"] [W/m2]")+"     "+("Preciznost ozračenja: ["+str(IRRADIANCE_LOWER_BOUND_FACTOR)+","+str(IRRADIANCE_UPPER_BOUND_FACTOR)+"]")
        _xLabel="Ozračenje FN modula,E [W/m2]"
        _yLabel="Učinkovitost FN modula,η"

        for j in _x:
            _xx.append(j)
        for j in _y:
            _yy.append(j)
            
    global AXIS_VALUE_PRECISION
    AXIS_VALUE_PRECISION=4
    
    return _xx,_yy,_xLabel,_yLabel,_details


def plot_mode_4(month, day, irradiance_lower_bound_factor, irradiance_upper_bound_factor, irradiance_values_start, irradiance_values_end, irradiance_values_step,temperature_mode,temperature_val, temperature_lower_bound_factor, temperature_upper_bound_factor):

    _xx=[]
    _yy=[]
    _x=[]
    _y=[]
    _partialModeText=""
    _modeTextString=""
    _details=""
    _xLabel=""
    _yLabel=""
    
    for irradiance_value in range(int(irradiance_values_start),int(irradiance_values_end),int(irradiance_values_step)):
        _x=get_data_with_set_irradiance(get_data_path(mm,dd,2), float(irradiance_value),irradiance_lower_bound_factor, irradiance_upper_bound_factor)
        _y,_partialModeText=get_data_with_set_temperature(get_data_path(mm,dd,1),temperature_mode, temperature_val, temperature_lower_bound_factor, temperature_upper_bound_factor)
        _details=("Ozračenje: ["+str(irradiance_values_start)+"-"+str(irradiance_values_end)+"] [W/m2]")+"     "+("Preciznost ozračenja: ["+str(IRRADIANCE_LOWER_BOUND_FACTOR)+","+str(IRRADIANCE_UPPER_BOUND_FACTOR)+"]")
        _xLabel="Ozračenje FN modula,E [W/m2]"
        _yLabel=_partialModeText.split("=")[1]

        for j in _x:
            _xx.append(j)
        for j in _y:
            _yy.append(j)
            
    global AXIS_VALUE_PRECISION
    AXIS_VALUE_PRECISION=0
    return _xx,_yy,_xLabel,_yLabel,_details

def get_data_for_two_parameter_dependency(month,day):

    _x=[]
    _y=[]
    _z=[]

    _xx=[]
    _yy=[]
    _zz=[]

    _partialModeText=""
    _modeText=""
    _xLabel=""
    _yLabel=""
    _zLabel=""

    for irradiance_value in range(int(SET_IRRADIANCE_VALUES[0]),int(SET_IRRADIANCE_VALUES[1]),int(SET_IRRADIANCE_VALUES[2])):
        _x,_partialModeText=get_data_with_set_temperature(get_data_path(month, day,1), TEMPERATURE_MODE, SET_TEMPERATURE, TEMPERATURE_LOWER_BOUND_FACTOR, TEMPERATURE_UPPER_BOUND_FACTOR)
        _y=get_data_with_set_irradiance(get_data_path(month, day,2), irradiance_value,IRRADIANCE_LOWER_BOUND_FACTOR, IRRADIANCE_UPPER_BOUND_FACTOR)
        
        if(PLOT_MODE==0):
            _z,_modeText=get_electrical_data_avg(get_data_path(month, day,0), _y, ELECTRICAL_PARAM_MODE)
        elif(PLOT_MODE==1):
            _z,_modeText=get_efficiency_irradiance(month, day,_y)
            
        else:
            _z.append(0)
            print("Nije dostupan mod!")

        for i in _x:
            _xx.append(i)
        for i in _y:
            _yy.append(i)
        for i in _z:
            _zz.append(i)
    
    _xLabel=_partialModeText.split("=")[1]
    _yLabel="Ozračenje FN modula,E [W/m2]"
    _zLabel=_modeText
        
    return _xx,_yy,_zz,_xLabel,_yLabel,_zLabel

def seek_maximum_efficiency(dataOfX,dataOfY,dataOfZ):
    _index_of_max_val=dataOfZ.index(max(dataOfZ))
    _result_temperature=dataOfX[_index_of_max_val]
    _result_irradiance=dataOfY[_index_of_max_val]
    _result_efficiency=dataOfZ[_index_of_max_val]

    return _result_temperature,_result_irradiance,_result_efficiency

def seek_index_of_maximum_efficiency(dataOfX,dataOfY,dataOfZ,dataOfTime):
    _index_of_max_val=-1

##    if((len(dataOfX)==len(dataOfY)) and (len(dataOfY)==len(dataOfZ)) and len((dataOfZ)==len(dataOfTime))):
    _index_of_max_val=dataOfZ.index(max(dataOfZ))
        
    return _index_of_max_val

def adjust_value_axis_processor(ax,max_dataOfY,typeOfAxis):
    
     if((max_dataOfY)):
        if(ADJUSTABLE_Y_AXIS==1 and SHOW_ONLY_LSR_LINE==False):
            _factor=0
            _values=[]
            _values.append(0)
            _values_max=float(max(max_dataOfY))
    
            for i in range(1,12):
                _factor+=0.1
                _values.append(round(_factor*_values_max,AXIS_VALUE_PRECISION))

            if(typeOfAxis==1):
                ax.set_yticks(_values[::1])
                ax.set_yticklabels(_values[::1])
                ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.'+str(AXIS_VALUE_PRECISION)+'f'))
            elif(typeOfAxis==2):
                ax.set_zticks(_values[::1])
                ax.set_zticklabels(_values[::1])
                ax.zaxis.set_major_formatter(plt.FormatStrFormatter('%.'+str(AXIS_VALUE_PRECISION)+'f'))

def get_dependent_histograms(fig,ax,hist_counts,hist_bins,mode,temperature_array,irradiance_array,efficiency_array,date_array,start_date,end_date):
    _max_counts_hist=max(hist_counts)
    _max_counts_hist_index=[]

    for i in range(len(hist_counts)):
        if (_max_counts_hist==hist_counts[i]):
            _max_counts_hist_index.append(int(i))

    _max_bins_hist=[]

    for i in _max_counts_hist_index:
        _max_bins_hist.append(hist_bins[i])

    _bin_width=hist_bins[1]-hist_bins[0]

    _bin_values_hist=[]

    for i in _max_counts_hist_index:
        _bin_values_hist.append(hist_bins[i])

    _temperature_values=[]
    _irradiance_values=[]
    _date_values=[]
    _efficiency_values=[]

    _reference_array=[]

    if(mode==0):
        _reference_array=temperature_array        
    elif(mode==1):
        _reference_array=irradiance_array        
            
    for i in range(len(_reference_array)):
        for j in _bin_values_hist:
            if(_reference_array[i]>=j and _reference_array[i]<=j+_bin_width):
                _temperature_values.append(temperature_array[i])
                _irradiance_values.append(irradiance_array[i])
                _date_values.append(date_array[i])
                _efficiency_values.append(efficiency_array[i])

    if(mode==0):
        ax.hist(_irradiance_values,color="green",ec="black",lw=1)
        ax.set_xlabel("Ozračenje [W/m2]",fontsize=LABEL_FONT_SIZE)

        if(end_date==start_date):
            plt.title('Mjerenja na dan '+start_date+'\n'+'Ozračenja za postizanje najveće '+DICT_PLOTLABEL[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))]+ ' u zavisnosti o utvrđenoj temperaturi',fontsize=LABEL_FONT_SIZE)
        else: plt.title('Mjerenja: '+start_date+' - '+end_date+'\n'+'Ozračenja za postizanje najveće '+DICT_PLOTLABEL[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))]+ ' u zavisnosti o utvrđenoj temperaturi',fontsize=LABEL_FONT_SIZE)
    elif(mode==1):
        ax.hist(_temperature_values,color="orange",ec="black",lw=1)
        ax.set_xlabel("Temperatura [oC]",fontsize=LABEL_FONT_SIZE)

        if(end_date==start_date):
            plt.title('Mjerenja na dan '+start_date+'\n'+'Temperatura za postizanje najveće '+DICT_PLOTLABEL[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))]+ ' u zavisnosti o utvrđenom ozračenju',fontsize=LABEL_FONT_SIZE)
        else: plt.title('Mjerenja: '+start_date+' - '+end_date+'\n'+'Temperatura za postizanje najveće '+DICT_PLOTLABEL[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))]+ ' u zavisnosti o utvrđenom ozračenju',fontsize=LABEL_FONT_SIZE)

    ax.set_ylabel("Frekvencija",fontsize=LABEL_FONT_SIZE)
                        
    if(SHOW_DATA_FOR_EXPORT==True):
        print('\nNajčešće temperature i preko njih ozračenja koje su se pojavljivale u skupu temperatura maksimalno postignutih učinkovitosti su:')
        print("DATUM\tTEMPERATURA [oC]\tOZRAČENJE [W/m2]\t"+DICT_PLOTLABEL[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].upper())
        for i in range(len(_temperature_values)):
            print(str(_date_values[i])+"\t"+str(_temperature_values[i])+"\t"+str(_irradiance_values[i])+"\t"+str(_efficiency_values[i]))
                
def additional_mode_1(ax):

    _dataOfX_2=[]
    _dataOfY_2=[]
    _dataOfZ_2=[]

    _partialDataOfX_2=[]
    _partialDataOfY_2=[]
    _partialDataOfZ_2=[]

    _xLabel_2=""
    _yLabel_2=""
    _zLabel_2=""

    _startDate=""
    _endDate=""
    _legendString=""

    _xx=[]
    _yy=[]
    _zz=[]

    _zz_max=[]

    _temperature_of_max_efficiency=[]
    _irradiance_of_max_efficiency=[]
    _max_efficiency=[]
    _date_of_max_efficiency=[]

    _temperature_of_max_efficiency_2=[]
    _irradiance_of_max_efficiency_2=[]
    _max_efficiency_2=[]
    _date_of_max_efficiency_2=[]

    _temperature_of_max_efficiency_monthly=[]
    _irradiance_of_max_efficiency_monthly=[]
    _max_efficiency_monthly=[]
    _date_of_max_efficiency_monthly=[]
    
    _efficiency_report=""

    for mm in range(int(MONTHS[0]),int(MONTHS[1])+1):    
        for dd in range(int(DAYS[0]),int(DAYS[1])+1):
        
            try:

                if(not(date_exists(mm,dd))):
                    _endDate=str(dd-1)+"."+str(mm)+"."
                    break
                else: _endDate=str(dd)+"."+str(mm)+"."

                _startDate=str(DAYS[0])+"."+str(mm)+"."

                _xx,_yy,_zz,_xLabel_2,_yLabel_2,_zLabel_2=get_data_for_two_parameter_dependency(mm,dd)

                _partialDataOfX_2,_partialDataOfY_2,_partialDataOfZ_2=match_XYZ_plot_data(_xx,_yy,_zz)

                for i in range(len(_partialDataOfX_2)):  
                    _dataOfX_2.append(_partialDataOfX_2[i])
                    _dataOfY_2.append(_partialDataOfY_2[i])
                    _dataOfZ_2.append(_partialDataOfZ_2[i])

                if(REMOVE_FAULT_DATA==True):
##                    REMOVE_FAULT_DATA==True and ADDITIONAL_MODE==3
                    remove_fault_data_XYZ(_partialDataOfX_2,_partialDataOfY_2,_partialDataOfZ_2)

                if(ADDITIONAL_MODE==3):

                    _temperature_of_max_efficiency_temp,_irradiance_of_max_efficiency_temp,_max_efficiency_temp=seek_maximum_efficiency(_partialDataOfX_2,_partialDataOfY_2,_partialDataOfZ_2)

                    if((_temperature_of_max_efficiency_temp)and(_irradiance_of_max_efficiency_temp)and(_max_efficiency_temp)):
                        _temperature_of_max_efficiency.append(float(_temperature_of_max_efficiency_temp))
                        _irradiance_of_max_efficiency.append(float(_irradiance_of_max_efficiency_temp))
                        _max_efficiency.append(float(_max_efficiency_temp))
                        _date_of_max_efficiency.append(str(dd)+"."+str(mm)+".")

                        _temperature_of_max_efficiency_2.append(float(_temperature_of_max_efficiency_temp))
                        _irradiance_of_max_efficiency_2.append(float(_irradiance_of_max_efficiency_temp))
                        _max_efficiency_2.append(float(_max_efficiency_temp))
                        _date_of_max_efficiency_2.append(str(dd)+"."+str(mm)+".")

                        _efficiency_report+=str(dd)+"."+str(mm)+".: Temperatura [oC]: "+str(_temperature_of_max_efficiency_temp)+" Ozračenje [W/m2]: "+str(_irradiance_of_max_efficiency_temp)+DICT_PLOTLABEL[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))]+str(_max_efficiency_temp)+"\n"

                if(MULTIPLE_DAYS_GRAPH==True and MULTIPLE_MONTHS_GRAPH==False and (_dataOfX_2)and (_dataOfY_2)and (_dataOfZ_2)):

                    if(REMOVE_FAULT_DATA==True):
                        remove_fault_data_XYZ(_dataOfX_2,_dataOfY_2,_dataOfZ_2)
                    
                    ax.scatter(_dataOfX_2,_dataOfY_2,_dataOfZ_2,label=str(dd)+"."+str(mm)+".")

                    if(ADJUSTABLE_Y_AXIS==1 and (_dataOfZ_2)):
                        _zz_max.append(float(max(_dataOfZ_2)))

                    ax.legend(loc=2)

                    _dataOfX_2=[]
                    _dataOfY_2=[]
                    _dataOfZ_2=[]
                pass

            except Exception as e:
                if(SHOW_LOG==1): print(e)

        if(_endDate==_startDate):
            _legendString=_startDate
        else: _legendString=_startDate+' - '+_endDate

        if(MULTIPLE_MONTHS_GRAPH==True and MULTIPLE_DAYS_GRAPH==False and (_dataOfX_2)and (_dataOfY_2)and (_dataOfZ_2) ):

            if(ADDITIONAL_MODE==3):
                _index_max_efficiency=seek_index_of_maximum_efficiency(_temperature_of_max_efficiency_2,_irradiance_of_max_efficiency_2,_max_efficiency_2,_date_of_max_efficiency_2)

                _temperature_of_max_efficiency_monthly.append(float(_temperature_of_max_efficiency_2[_index_max_efficiency]))
                _irradiance_of_max_efficiency_monthly.append(float(_irradiance_of_max_efficiency_2[_index_max_efficiency]))
                _max_efficiency_monthly.append(float(_max_efficiency_2[_index_max_efficiency]))
                _date_of_max_efficiency_monthly.append(str(_date_of_max_efficiency_2[_index_max_efficiency]))
                #DICT_DATE_MONTH[str(_date_of_max_efficiency_2[_index_max_efficiency]).split(".")[1]]

                _temperature_of_max_efficiency_2=[]
                _irradiance_of_max_efficiency_2=[]
                _max_efficiency_2=[]
                _date_of_max_efficiency_2=[]

            if(REMOVE_FAULT_DATA==True):
                remove_fault_data_XYZ(_dataOfX_2,_dataOfY_2,_dataOfZ_2)
            
            ax.scatter(_dataOfX_2,_dataOfY_2,_dataOfZ_2,label=_legendString)

            if(ADJUSTABLE_Y_AXIS==1 and (_dataOfZ_2)):
                _zz_max.append(float(max(_dataOfZ_2)))

            ax.legend(loc=2)

            _dataOfX_2=[]
            _dataOfY_2=[]
            _dataOfZ_2=[]

    if(MULTIPLE_MONTHS_GRAPH==False and MULTIPLE_DAYS_GRAPH==False and (_dataOfX_2)and (_dataOfY_2)and (_dataOfZ_2)):
        ##(MULTIPLE_MONTHS_GRAPH==False and MULTIPLE_DAYS_GRAPH==False and ADDITIONAL_MODE==1) or ADDITIONAL_MODE==3
        if(REMOVE_FAULT_DATA==True):
            remove_fault_data_XYZ(_dataOfX_2,_dataOfY_2,_dataOfZ_2)
            
        ax.scatter(_dataOfX_2,_dataOfY_2,_dataOfZ_2)

        if(ADJUSTABLE_Y_AXIS==1 and (_dataOfZ_2)):
            _zz_max.append(float(max(_dataOfZ_2)))

    if(ADJUSTABLE_Y_AXIS==1):
        adjust_value_axis_processor(ax,_zz_max,2)    

    _startDate=str(DAYS[0])+"."+str(MONTHS[0])+"."

    if(_endDate==_startDate):
        plt.title('Mjerenja na dan '+_startDate+'\n',fontsize=LABEL_FONT_SIZE)
    else: plt.title('Mjerenja: '+_startDate+' - '+_endDate+'\n',fontsize=LABEL_FONT_SIZE)

    ax.set_xlabel("\n\n\n"+_xLabel_2,fontsize=LABEL_FONT_SIZE)
    ax.set_ylabel("\n\n\n"+_yLabel_2,fontsize=LABEL_FONT_SIZE)
    ax.set_zlabel("\n\n\n"+_zLabel_2,fontsize=LABEL_FONT_SIZE)

    ax.tick_params(axis='both', which='major', labelsize=LABEL_FONT_SIZE)

    if(SHOW_STATISTICS==True):
        print("Temperatura i ozračenje za najbolju "+DICT_PLOTLABEL[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))]+": ")
        print(_efficiency_report)

    if(ADDITIONAL_MODE==3):
            
        _fig,_ax=plt.subplots()

        show_fault_data_XYZDate(_temperature_of_max_efficiency,_irradiance_of_max_efficiency,_max_efficiency,_date_of_max_efficiency)
        additional_mode_3(_ax,_temperature_of_max_efficiency,_irradiance_of_max_efficiency,_max_efficiency,_date_of_max_efficiency,_endDate,_zLabel_2)

        if(MULTIPLE_DAYS_GRAPH==False and MULTIPLE_MONTHS_GRAPH==False and (_temperature_of_max_efficiency) and (_irradiance_of_max_efficiency)):

            _fig2,_ax2=plt.subplots()
            _counts,_bins,_bars=_ax2.hist(_temperature_of_max_efficiency,color="orange",ec="black",lw=1)

            if(_endDate==_startDate):
                plt.title('Mjerenja na dan '+_startDate+'\n'+'Temperature za postizanje najveće '+DICT_PLOTLABEL[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))],fontsize=LABEL_FONT_SIZE)
            else: plt.title('Mjerenja: '+_startDate+' - '+_endDate+'\n'+'Temperature za postizanje najveće '+DICT_PLOTLABEL[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))],fontsize=LABEL_FONT_SIZE)

            _ax2.set_xlabel("Temperatura [oC]",fontsize=LABEL_FONT_SIZE)
            _ax2.set_ylabel("Frekvencija",fontsize=LABEL_FONT_SIZE)

            _fig3,_ax3=plt.subplots()

            get_dependent_histograms(_fig3,_ax3,_counts,_bins,0,_temperature_of_max_efficiency,_irradiance_of_max_efficiency,_max_efficiency,_date_of_max_efficiency,_startDate,_endDate)
        
            _fig4,_ax4=plt.subplots()
            _counts,_bins,_bars=_ax4.hist(_irradiance_of_max_efficiency,color="green",ec="black",lw=1)
        
            if(_endDate==_startDate):
                plt.title('Mjerenja na dan '+_startDate+'\n'+'Ozračenja za postizanje najveće '+DICT_PLOTLABEL[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))],fontsize=LABEL_FONT_SIZE)
            else: plt.title('Mjerenja: '+_startDate+' - '+_endDate+'\n'+'Ozračenja za postizanje najveće '+DICT_PLOTLABEL[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))],fontsize=LABEL_FONT_SIZE)
        
            _ax4.set_xlabel("Ozračenje [W/m2]",fontsize=LABEL_FONT_SIZE)
            _ax4.set_ylabel("Frekvencija",fontsize=LABEL_FONT_SIZE)

            _fig5,_ax5=plt.subplots()
            
            get_dependent_histograms(_fig5,_ax5,_counts,_bins,1,_temperature_of_max_efficiency,_irradiance_of_max_efficiency,_max_efficiency,_date_of_max_efficiency,_startDate,_endDate)

        if(MULTIPLE_MONTHS_GRAPH==True and MULTIPLE_DAYS_GRAPH==False and (_date_of_max_efficiency_monthly)and (_max_efficiency_monthly)):
        
            _fig4,_ax4=plt.subplots()

            _x=[dt.datetime.strptime(date.split(".")[1], '%m') for date in _date_of_max_efficiency_monthly]

            plt.gcf().autofmt_xdate()
        
            plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%B'))
            plt.gca().xaxis.set_major_locator(mdates.MonthLocator())            
            
            _ax4.scatter(_x,_max_efficiency_monthly)
            
            if(_endDate==_startDate):
                plt.title('Mjerenja na dan '+_startDate+'\n'+'Najveća '+DICT_PLOTLABEL[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))]+' za promatrani mjesec',fontsize=LABEL_FONT_SIZE)
            else: plt.title('Mjerenja: '+_startDate+' - '+_endDate+'\n'+'Najveća '+DICT_PLOTLABEL[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))]+' za promatrani mjesec',fontsize=LABEL_FONT_SIZE)

            _ax4.set_xlabel("Mjesec",fontsize=LABEL_FONT_SIZE)
            _ax4.set_ylabel(_zLabel_2,fontsize=LABEL_FONT_SIZE)

            if(ADJUSTABLE_Y_AXIS==1):
                adjust_value_axis_processor(_ax4,_max_efficiency_monthly,1)

            if(SHOW_DATA_FOR_EXPORT==True):
                print("\nPodaci za izradu tablice  [Najveća "+DICT_PLOTLABEL[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))]+" za određeni mjesec]:")
                print("MJESEC\tDATUM\tTEMPERATURA [oC]\tOZRAČENJE [W/m2]\t"+DICT_PLOTLABEL[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].upper())
                for i in range(len(_date_of_max_efficiency_monthly)):
                    print(DICT_DATE_MONTH[str(_date_of_max_efficiency_monthly[i]).split(".")[1]]+"\t"+str(_date_of_max_efficiency_monthly[i])+"\t"+str(_temperature_of_max_efficiency_monthly[i])+"\t"+str(_irradiance_of_max_efficiency_monthly[i])+"\t"+str(_max_efficiency_monthly[i]))

            



def additional_mode_2(ax):

    _x=[]
    _y=[]

    _xx=[]
    _yy=[]

    _yy_max=[]

    _coeffs=[]

    _partialModeText=""
    _modeText=""
    _xLabel=""
    _yLabel=""

    _startDate=""
    _endDate=""

    _LSR_details=""
    _LSR_details_selected=""

    _best_efficiency_slope_original={}
    _best_efficiency_slope={}
    _best_efficiency_offset={}

    _best_efficiency_slope_original_selected={}
    _best_efficiency_slope_selected={}
    _best_efficiency_offset_selected={}


    _startDate=str(DAYS[0])+"."+str(MONTHS[0])+"."

    for irradiance_value in range(int(SET_IRRADIANCE_VALUES[0]),int(SET_IRRADIANCE_VALUES[1]),int(SET_IRRADIANCE_VALUES[2])):

        for mm in range(int(MONTHS[0]),int(MONTHS[1])+1):    
            for dd in range(int(DAYS[0]),int(DAYS[1])+1):
    
                if(not(date_exists(mm,dd))):
                    _endDate=str(dd-1)+"."+str(mm)+"."
                    break
                else: _endDate=str(dd)+"."+str(mm)+"."

                _x,_partialModeText=get_data_with_set_temperature(get_data_path(mm, dd,1), TEMPERATURE_MODE, SET_TEMPERATURE, TEMPERATURE_LOWER_BOUND_FACTOR, TEMPERATURE_UPPER_BOUND_FACTOR)
                _z=get_data_with_set_irradiance(get_data_path(mm, dd,2), irradiance_value,IRRADIANCE_LOWER_BOUND_FACTOR, IRRADIANCE_UPPER_BOUND_FACTOR)

                if(PLOT_MODE==0):
                    _y,_modeText=get_electrical_data_avg(get_data_path(mm, dd,0), _z, ELECTRICAL_PARAM_MODE)

                elif(PLOT_MODE==1):
                    _y,_modeText=get_efficiency_irradiance(mm,dd,_z)

                else:
                    _y.append(0)
                    print("Nije dostupan mod!")

                _x,_y=match_XY_plot_data(_x,_y)

                for i in _x:
                    _xx.append(i)
                for i in _y:
                    _yy.append(i)

        if(REMOVE_FAULT_DATA==True):
            remove_fault_data(_xx,_yy)

        if(SHOW_ONLY_LSR_LINE==False):
            ax.scatter(_xx,_yy,label="E="+str(irradiance_value)+" W/m2")
            _coeffs=fit_LSR(_xx,_yy,SHOW_LSR_LINE,SHOW_ONLY_LSR_LINE, 1,str("E="+str(irradiance_value)+" W/m2"))
        else:
            _coeffs=fit_LSR(_xx,_yy,1,SHOW_ONLY_LSR_LINE,1,str("E="+str(irradiance_value)+" W/m2"))

        _LSR_details+="E="+str(irradiance_value)+" W/m2 : "+DICT_PLOTMODE[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].split(";")[1]+"=["+str(_coeffs[0])+"]*"+DICT_PLOTMODE[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].split(";")[0]+"+["+str(_coeffs[1])+"];  "+"r-koeficijent: "+str(_coeffs[2])+"\n"                

        if((_coeffs[0]!=0.0)  and (_coeffs[1]!=0.0)):
            _best_efficiency_slope_original[irradiance_value]=_coeffs[0]
            _best_efficiency_slope[irradiance_value]=abs(_coeffs[0])
            _best_efficiency_offset[irradiance_value]=_coeffs[1]

            _best_efficiency_slope=dict(sorted(_best_efficiency_slope.items(), key=lambda item: item[1]))
            _best_efficiency_offset=dict(sorted(_best_efficiency_offset.items(), key=lambda item: item[1]))
        
        
        if(_coeffs[2]<=-0.7 or _coeffs[2]>=0.7 and ((_coeffs[0]!=0.0)  and (_coeffs[1]!=0.0))):
            _LSR_details_selected+="E="+str(irradiance_value)+" W/m2 : "+DICT_PLOTMODE[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].split(";")[1]+"=["+str(_coeffs[0])+"]*"+DICT_PLOTMODE[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].split(";")[0]+"+["+str(_coeffs[1])+"];  "+"r-koeficijent: "+str(_coeffs[2])+"\n"                

            _best_efficiency_slope_original_selected[irradiance_value]=_coeffs[0]
            _best_efficiency_slope_selected[irradiance_value]=abs(_coeffs[0])
            _best_efficiency_offset_selected[irradiance_value]=_coeffs[1]

            _best_efficiency_slope_selected=dict(sorted(_best_efficiency_slope_selected.items(), key=lambda item: item[1]))
            _best_efficiency_offset_selected=dict(sorted(_best_efficiency_offset_selected.items(), key=lambda item: item[1]))
            
        if(irradiance_value==200):
            print('Ekstrapolacija vrijednosti za E = 200 W/m2 pri 25 oC: '+str(_coeffs[0]*25.0+_coeffs[1]))
        elif(irradiance_value==1000):
            print('Ekstrapolacija vrijednosti za E = 1000 W/m2 pri 25 oC: '+str(_coeffs[0]*25.0+_coeffs[1]))

        if(ADJUSTABLE_Y_AXIS==1 and (_yy)):
            _yy_max.append(float(max(_yy)))

        _xx=[]
        _yy=[]

        

    _xLabel=_partialModeText.split("=")[1]
    _yLabel=_modeText

    plt.xlabel(_xLabel,fontsize=LABEL_FONT_SIZE)
    plt.ylabel(_yLabel,fontsize=LABEL_FONT_SIZE)

    ax.tick_params(axis='both', which='major', labelsize=LABEL_FONT_SIZE)

    if(_endDate==_startDate):
        plt.title('Mjerenja na dan '+_endDate+'\n',fontsize=LABEL_FONT_SIZE)
    else: plt.title('Mjerenja: '+_startDate+"-"+_endDate+'\n',fontsize=LABEL_FONT_SIZE)

    plt.legend(loc='upper right')

    print("\nJednadžba LSR:\n"+_LSR_details)
    print("Najpovoljniji temperaturni koeficijent FN modula je na E="+str(next(iter(_best_efficiency_slope)))+" W/m2: "+str(_best_efficiency_slope_original[next(iter(_best_efficiency_slope))])+" 1/K")
    print("Najpovoljnija "+DICT_PLOTLABEL[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))]+" FN modula na E="+str(next(reversed(_best_efficiency_offset.keys()))  )+" W/m2: "+str(_best_efficiency_offset[next(reversed(_best_efficiency_offset.keys()))]))
    print("Najnepovoljniji temperaturni koeficijent FN modula na E="+str(next(reversed(_best_efficiency_slope.keys()))  )+" W/m2: "+str(_best_efficiency_slope_original[next(reversed(_best_efficiency_slope.keys()))])+" 1/K")
    print("Najnepovoljnija "+DICT_PLOTLABEL[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))]+" FN modula je na E="+str(next(iter(_best_efficiency_offset)))+" W/m2: "+str(_best_efficiency_offset[next(iter(_best_efficiency_offset))])+"\n")


    print("\nJednadžba LSR koje su dovoljno dobro fit-ane na mjerene podatke:\n"+_LSR_details_selected)
    print("Najpovoljniji temperaturni koeficijent FN modula je na E="+str(next(iter(_best_efficiency_slope_selected)))+" W/m2: "+str(_best_efficiency_slope_original_selected[next(iter(_best_efficiency_slope_selected))])+" 1/K")
    print("Najpovoljnija "+DICT_PLOTLABEL[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))]+" FN modula na E="+str(next(reversed(_best_efficiency_offset_selected.keys()))  )+" W/m2: "+str(_best_efficiency_offset_selected[next(reversed(_best_efficiency_offset_selected.keys()))]))
    print("Najnepovoljniji temperaturni koeficijent FN modula na E="+str(next(reversed(_best_efficiency_slope_selected.keys()))  )+" W/m2: "+str(_best_efficiency_slope_original_selected[next(reversed(_best_efficiency_slope_selected.keys()))])+" 1/K")
    print("Najnepovoljnija "+DICT_PLOTLABEL[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))]+" FN modula je na E="+str(next(iter(_best_efficiency_offset_selected)))+" W/m2: "+str(_best_efficiency_offset_selected[next(iter(_best_efficiency_offset_selected))])+"\n")

    if(ADJUSTABLE_Y_AXIS==1):
        adjust_value_axis_processor(ax,_yy_max,1)

    
def additional_mode_3(ax,temperature_of_max_efficiency,irradiance_of_max_efficiency,max_efficiency,date_of_max_efficiency,_end_date, label_mode):

    if(ADJUSTABLE_X_AXIS == True):
        _x=[dt.datetime.strptime(date, '%d.%m.').date() for date in date_of_max_efficiency]
        
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d.%m.'))
        plt.gca().xaxis.set_major_locator(mdates.DayLocator())
        _interval_mode=1
        
        if(len(_x)>31 or (MULTIPLE_DAYS_GRAPH==False and MULTIPLE_MONTHS_GRAPH==True)):
            _interval=float(len(_x)/31.0)
            _interval_mode=TIME_SCALE*math.ceil(_interval)

        plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=_interval_mode))             
        ax.scatter(_x,max_efficiency)
        plt.gcf().autofmt_xdate()

    else:
        ax.scatter(date_of_max_efficiency,max_efficiency)

    global AXIS_VALUE_PRECISION
    AXIS_VALUE_PRECISION=4

    if(ADJUSTABLE_Y_AXIS==1):
        adjust_value_axis_processor(ax,max_efficiency,1)

    _startDate=str(DAYS[0])+"."+str(MONTHS[0])+"."

    if(_end_date==_startDate):
        plt.title('Mjerenja na dan '+_startDate+'\n',fontsize=LABEL_FONT_SIZE)
    else: plt.title('Mjerenja: '+_startDate+' - '+_end_date+'\n',fontsize=LABEL_FONT_SIZE)

    ax.set_xlabel("Datum",fontsize=LABEL_FONT_SIZE)
    ax.set_ylabel(label_mode,fontsize=LABEL_FONT_SIZE)

    if(SHOW_DATA_FOR_EXPORT==True):
        print("\nPodaci za izradu tablice [Najveća "+DICT_PLOTLABEL[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))]+" za određeni dan]:")
        print("DATUM\tTEMPERATURA [oC]\tOZRAČENJE [W/m2]\t"+DICT_PLOTLABEL[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].upper())
        for i in range(len(date_of_max_efficiency)):
            print(str(date_of_max_efficiency[i])+"\t"+str(temperature_of_max_efficiency[i])+"\t"+str(irradiance_of_max_efficiency[i])+"\t"+str(max_efficiency[i]))

def additional_mode_4(fig,ax):

    _irradiance=[]

    _values=[]
    _x=[]
    _y=[]
    
    _partialModeText=""

    _startDate=""
    _endDate=""
    _legendString=""

    global AXIS_VALUE_PRECISION

    for mm in range(int(MONTHS[0]),int(MONTHS[1])+1):
        for dd in range(int(DAYS[0]),int(DAYS[1])+1):
        
            try:

                if(not(date_exists(mm,dd))):
                    _endDate=str(dd-1)+"."+str(mm)+"."
                    break
                else: _endDate=str(dd)+"."+str(mm)+"."

                _startDate=str(DAYS[0])+"."+str(mm)+"."

                _irradiance=get_data_with_set_irradiance(get_data_path(mm,dd,2), -1,IRRADIANCE_LOWER_BOUND_FACTOR, IRRADIANCE_UPPER_BOUND_FACTOR)

                if(TIME_MODE==0):
                    _values=_irradiance
                    _partialModeText="Ozračenje,E [W/m2]"
                    AXIS_VALUE_PRECISION=2
                elif(TIME_MODE==1):
                    _values,_partialModeText=get_data_with_set_temperature(get_data_path(mm,dd,1),TEMPERATURE_MODE, -1, TEMPERATURE_LOWER_BOUND_FACTOR, TEMPERATURE_UPPER_BOUND_FACTOR)
                    _partialModeText=_partialModeText.split("=")[1]
                    AXIS_VALUE_PRECISION=2
                elif(TIME_MODE==2):
                    _values,_partialModeText=get_electrical_data_avg(get_data_path(mm,dd,0), _irradiance, ELECTRICAL_PARAM_MODE)
                    AXIS_VALUE_PRECISION=2
                elif(TIME_MODE==3):
                    _values,_partialModeText=get_efficiency_irradiance(mm,dd,_irradiance)
                    AXIS_VALUE_PRECISION=4
        
                for i in range(1,len(_values),2):
                    _x.append((_values[i])[11:16])

                for i in range(0,len(_values),2):
                    _y.append(float(_values[i]))
                        
                if(MULTIPLE_DAYS_GRAPH==True and MULTIPLE_MONTHS_GRAPH==False and (_x) and (_y)):
                    _xx=[dt.datetime.strptime(time, '%H:%M') for time in _x]
                    ax.plot(_xx,_y,label=str(dd)+"."+str(mm)+".")

                    plt.legend(loc='upper right')

                    _x=[]
                    _xx=[]
                    _y=[]

                pass

            except Exception as e:
                print(e)

    ##endDate=str(dd)+"."+str(mm)+"."
               
        if(_endDate==_startDate):
            _legendString=_startDate
        else: _legendString=_startDate+' - '+_endDate

        if(MULTIPLE_MONTHS_GRAPH==True and MULTIPLE_DAYS_GRAPH==False and (_x) and (_y)):
            _xx=[dt.datetime.strptime(time, '%H:%M') for time in _x]
            ax.plot(_xx,_y,label=_legendString)

            plt.legend(loc='upper right')

            _x=[]
            _xx=[]
            _y=[]    
            
    if(MULTIPLE_MONTHS_GRAPH==False and MULTIPLE_DAYS_GRAPH==False and (_x) and (_y)):
        _xx=[dt.datetime.strptime(time, '%H:%M') for time in _x]
        ax.plot(_xx,_y)

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.gca().xaxis.set_major_locator(mdates.HourLocator())        
    plt.gca().xaxis.set_major_locator(mdates.HourLocator(interval=2))
    plt.gcf().autofmt_xdate()

    ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.'+str(AXIS_VALUE_PRECISION)+'f'))

    if(_endDate==_startDate):
        plt.title('Mjerenja na dan '+_startDate+'\n',fontsize=LABEL_FONT_SIZE)
    else: plt.title('Mjerenja: '+_startDate+' - '+_endDate+'\n',fontsize=LABEL_FONT_SIZE)

    ax.set_xlabel('Vrijeme,t',fontsize=LABEL_FONT_SIZE)
    ax.set_ylabel(_partialModeText,fontsize=LABEL_FONT_SIZE)


##MAIN OF PROGRAM####################################################################################################

dataOfX=[]
dataOfY=[]

partialDataOfX=[]
partialDataOfY=[]

xx=[]
yy=[]

yy_max=[]

xLabel=""
yLabel=""

startDate=""
endDate=""

details=""
LSR_details=""
legendString=""

coeffs=[]
statistics_data=[]

##MAIN PLOT#####################################################################

locale.setlocale(locale.LC_ALL,'hr_HR')

fig,ax=plt.subplots()

print("\n"+str(DICT_TYPE_OF_PV_PANEL[TYPE_OF_PV_PANEL])+" MODUL\n")

for mm in range(int(MONTHS[0]),int(MONTHS[1])+1):
    for dd in range(int(DAYS[0]),int(DAYS[1])+1):
        
        try:

            if(not(date_exists(mm,dd))):
                endDate=str(dd-1)+"."+str(mm)+"."
                break
            else: endDate=str(dd)+"."+str(mm)+"."

            startDate=str(DAYS[0])+"."+str(mm)+"."
                        
            if PLOT_MODE==0:
                xx,yy,xLabel,yLabel,details=plot_mode_0(mm,dd, SET_IRRADIANCE, IRRADIANCE_LOWER_BOUND_FACTOR,IRRADIANCE_UPPER_BOUND_FACTOR, TEMPERATURE_MODE, SET_TEMPERATURE,TEMPERATURE_LOWER_BOUND_FACTOR,TEMPERATURE_UPPER_BOUND_FACTOR,ELECTRICAL_PARAM_MODE)           
            elif PLOT_MODE==1:
                xx,yy,xLabel,yLabel,details=plot_mode_1(mm,dd, SET_IRRADIANCE, IRRADIANCE_LOWER_BOUND_FACTOR,IRRADIANCE_UPPER_BOUND_FACTOR, TEMPERATURE_MODE, SET_TEMPERATURE,TEMPERATURE_LOWER_BOUND_FACTOR,TEMPERATURE_UPPER_BOUND_FACTOR)           
            elif PLOT_MODE==2:
                xx,yy,xLabel,yLabel,details=plot_mode_2(mm,dd, IRRADIANCE_LOWER_BOUND_FACTOR,IRRADIANCE_UPPER_BOUND_FACTOR, SET_IRRADIANCE_VALUES[0], SET_IRRADIANCE_VALUES[1], SET_IRRADIANCE_VALUES[2],ELECTRICAL_PARAM_MODE)
            elif PLOT_MODE==3:
                xx,yy,xLabel,yLabel,details=plot_mode_3(mm,dd, IRRADIANCE_LOWER_BOUND_FACTOR,IRRADIANCE_UPPER_BOUND_FACTOR, SET_IRRADIANCE_VALUES[0], SET_IRRADIANCE_VALUES[1], SET_IRRADIANCE_VALUES[2])                 
            elif PLOT_MODE==4:
                xx,yy,xLabel,yLabel,details=plot_mode_4(mm,dd, IRRADIANCE_LOWER_BOUND_FACTOR,IRRADIANCE_UPPER_BOUND_FACTOR, SET_IRRADIANCE_VALUES[0], SET_IRRADIANCE_VALUES[1], SET_IRRADIANCE_VALUES[2],TEMPERATURE_MODE, SET_TEMPERATURE,TEMPERATURE_LOWER_BOUND_FACTOR,TEMPERATURE_UPPER_BOUND_FACTOR)

            partialDataOfX,partialDataOfY=match_XY_plot_data(xx,yy)

            for i in range(len(partialDataOfX)):  
                dataOfX.append(partialDataOfX[i])
                dataOfY.append(partialDataOfY[i])

            if(len(partialDataOfX)>1 and len(partialDataOfY)>1):

                if(MULTIPLE_DAYS_GRAPH==True and MULTIPLE_MONTHS_GRAPH==False):

                    if(REMOVE_FAULT_DATA==True):
                        remove_fault_data(dataOfX,dataOfY)
    
                    if(SHOW_ONLY_LSR_LINE==False):
                        ax.scatter(dataOfX,dataOfY,label=str(dd)+"."+str(mm)+".")
                        coeffs=fit_LSR(dataOfX,dataOfY,SHOW_LSR_LINE,SHOW_ONLY_LSR_LINE,1,endDate)
                    else:
                        coeffs=fit_LSR(dataOfX,dataOfY,1,SHOW_ONLY_LSR_LINE,1,endDate)

                    LSR_details+=endDate+" : "+DICT_PLOTMODE[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].split(";")[1]+"=["+str(coeffs[0])+"]*"+DICT_PLOTMODE[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].split(";")[0]+"+["+str(coeffs[1])+"];  "+"r-koeficijent: "+str(coeffs[2])+"\n"             
                    
                    if(ADJUSTABLE_Y_AXIS==1 and (dataOfY)): #second statement in if checks if dataOfY is NOT EMPTY
                        yy_max.append(float(max(dataOfY)))
                    
                    if(SHOW_STATISTICS==True):
                        print("\n"+str(dd)+"."+str(mm)+".:")
                        statistics_data=get_statistics(dataOfX,dataOfY)
                        print("Temperaturni koeficijent d"+DICT_PLOTMODE[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].split(";")[1]+"/d"+DICT_PLOTMODE[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].split(";")[0]+": "+str(coeffs[0])+" 1/K \n")
                    
                    plt.legend(loc='upper right')

                    dataOfX=[]
                    dataOfY=[]

            pass

        except Exception as e:
            print(e)

    ##endDate=str(dd)+"."+str(mm)+"."
               
    if(endDate==startDate):
        legendString=startDate
    else: legendString=startDate+' - '+endDate

    if(len(dataOfX)>1 and len(dataOfY)>1):

        if(MULTIPLE_MONTHS_GRAPH==True and MULTIPLE_DAYS_GRAPH==False):

            if(REMOVE_FAULT_DATA==True):
                remove_fault_data(dataOfX,dataOfY)

            if(SHOW_ONLY_LSR_LINE==False):
                ax.scatter(dataOfX,dataOfY,label=legendString)
                coeffs=fit_LSR(dataOfX,dataOfY,SHOW_LSR_LINE,SHOW_ONLY_LSR_LINE,1,legendString)
            else:
                coeffs=fit_LSR(dataOfX,dataOfY,1,SHOW_ONLY_LSR_LINE,1,legendString)

            LSR_details+=str(DAYS[0])+"."+str(mm)+" - "+endDate+" : "+DICT_PLOTMODE[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].split(";")[1]+"=["+str(coeffs[0])+"]*"+DICT_PLOTMODE[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].split(";")[0]+"+["+str(coeffs[1])+"];  "+"r-koeficijent: "+str(coeffs[2])+"\n"                

            if(ADJUSTABLE_Y_AXIS==1 and (dataOfY)):
                yy_max.append(float(max(dataOfY)))

            if(SHOW_STATISTICS==True):
                print("\n"+str(DAYS[0])+"."+str(mm)+" - "+endDate+":")
                statistics_data=get_statistics(dataOfX,dataOfY)
                print("Temperaturni koeficijent d"+DICT_PLOTMODE[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].split(";")[1]+"/d"+DICT_PLOTMODE[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].split(";")[0]+": "+str(coeffs[0])+" 1/K \n")

                        
            plt.legend(loc='upper right')
            
            dataOfX=[]
            dataOfY=[]

if(len(dataOfX)>1 and len(dataOfY)>1):

    if(MULTIPLE_MONTHS_GRAPH==False and MULTIPLE_DAYS_GRAPH==False and (dataOfX) and (dataOfY)):

        if(REMOVE_FAULT_DATA==True):
            remove_fault_data(dataOfX,dataOfY)

        if(SHOW_ONLY_LSR_LINE==False):
            ax.scatter(dataOfX,dataOfY)
            coeffs=fit_LSR(dataOfX,dataOfY,SHOW_LSR_LINE,SHOW_ONLY_LSR_LINE,0,legendString)
        else:
            coeffs=fit_LSR(dataOfX,dataOfY,1,SHOW_ONLY_LSR_LINE,0,legendString)

        LSR_details+=str(DAYS[0])+"."+str(MONTHS[0])+". - "+endDate+" : "+DICT_PLOTMODE[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].split(";")[1]+"=["+str(coeffs[0])+"]*"+DICT_PLOTMODE[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].split(";")[0]+"+["+str(coeffs[1])+"];  "+"r-koeficijent: "+str(coeffs[2])+"\n"                
        yy_max.append(float(max(dataOfY)))

        if(SHOW_STATISTICS==True):
            print("\n"+str(DAYS[0])+"."+str(MONTHS[0])+" - "+endDate+":")
            statistics_data=get_statistics(dataOfX,dataOfY)
            print("Temperaturni koeficijent d"+DICT_PLOTMODE[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].split(";")[1]+"/d"+DICT_PLOTMODE[str(str(PLOT_MODE)+str(ELECTRICAL_PARAM_MODE))].split(";")[0]+": "+str(coeffs[0])+" 1/K \n")


startDate=str(DAYS[0])+"."+str(MONTHS[0])+"."
if(endDate==startDate):
    plt.title('Mjerenja na dan '+startDate+'\n',fontsize=LABEL_FONT_SIZE)
else: plt.title('Mjerenja: '+startDate+' - '+endDate+'\n',fontsize=LABEL_FONT_SIZE)

if(ADJUSTABLE_Y_AXIS==1):
    adjust_value_axis_processor(ax,yy_max,1)

if(SHOW_STATISTICS==True):
    print("\nJednadžba LSR:\n"+LSR_details)

plt.xlabel(xLabel+'\n\n'+details+'\n\n'+LSR_details,fontsize=LABEL_FONT_SIZE)
plt.ylabel(yLabel,fontsize=LABEL_FONT_SIZE)

ax.tick_params(axis='both', which='major', labelsize=LABEL_FONT_SIZE)


## ADDITIONAL PLOT ####################################################################################

if(ADDITIONAL_MODE==1 or ADDITIONAL_MODE==3 and (PLOT_MODE==0 or PLOT_MODE==1)):
    fig2=plt.figure()
    ax2 = fig2.add_subplot(projection='3d')
    additional_mode_1(ax2)
   
if(ADDITIONAL_MODE==2 and (PLOT_MODE==0 or PLOT_MODE==1)):    
    fig3,ax3=plt.subplots()    
    additional_mode_2(ax3)

if(ADDITIONAL_MODE==4):
    fig4,ax4=plt.subplots()
    additional_mode_4(fig4,ax4)


plt.show()
