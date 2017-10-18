from eppy import modeleditor
from eppy.modeleditor import IDF

iddfile = "C:/openstudio-2.2.0/EnergyPlus/Energy+.idd"
try:
  IDF.setiddname(iddfile)
except modeleditor.IDDAlreadySetError as e:
  pass

fname = "./Export_2_LargerTest.idf"
idf1 = IDF(fname)

# idf1.printidf()

# Directions that will be co-simulated
# Starting with South just to start but can expand to
# South, East, West, and Core
# Case sensitive, so match the case in the OpenStudio file
directions = ['South']
daylighting = True

instructions = """
For this script to run correctly: 
	1) Make sure that all components have a consistent name that is formatted by BITCoPT <direction> Component Name"
		a) Test components include, it is recommended to follow the same naming convention: 
			i) Fenestration:Detailed => BITCoPT <direction> Window
			ii) PlantComponent:TemperatureSource => BITCoPT <direction> Plant Component Temperature Source
			iii) Pump:Variable => BITCoPT <direction> Var Spd Pump
			iv) SetpointManager:Scheduled => BITCoPT <direction> Scheduled HW Temp
	2) If your file already includes an ExternalInterface component make sure it is for FunctionalMockupUnitExport
	3) If there are already generators in the base IDF file, change the conditional statement on line 200
	4) Make sure you have an individual lighting definition for each zone in OpenStudio so that they can be controlled independently with the daylighting schedule in the co-simulation
		a) Name each lighting schedule by the name of the zone and direction: South Lights
"""
print instructions

print 'For this script to run correctly: '
print

print "The IDF file you specified at "+str(fname)+" will be modified for co-simulation.\n"

print "The surface(s) being modified include: " + str(directions)+" \n"


print '_________________________________________________________________\n'

print "Adding the ExternalInterface:FunctionalMockupUnitExport component:"

# Create new object for the ExternalInterface needed for FMU export
if len(idf1.idfobjects['EXTERNALINTERFACE']) < 1:
	idf1.newidfobject('EXTERNALINTERFACE')
	idf1.idfobjects['EXTERNALINTERFACE'][0].Name_of_External_Interface = 'FunctionalMockupUnitExport'

print idf1.idfobjects['EXTERNALINTERFACE'][0]

print '_________________________________________________________________\n'

print "Adding the ExternalInterface:FunctionalMockupUnitExport:From:Variable component for the Outdoor Air Temperature:"

zones = idf1.idfobjects['ZONE']

# Add FMU Output object for Outdoor Air Temperature
idf1.newidfobject('EXTERNALINTERFACE:FUNCTIONALMOCKUPUNITEXPORT:FROM:VARIABLE')
fromVariableLast = idf1.idfobjects['EXTERNALINTERFACE:FUNCTIONALMOCKUPUNITEXPORT:FROM:VARIABLE'][-1]
fromVariableLast.OutputVariable_Index_Key_Name = zones[-1].Name
fromVariableLast.OutputVariable_Name = 'Zone Outdoor Air Drybulb Temperature'
fromVariableLast.FMU_Variable_Name = 'TOutEnv'

print fromVariableLast




# Setup windows objects for co-simulation
BITCoPTWindows = []
for name in directions:
	BITCoPTWindows.append('BITCoPT ' + str(name) + ' Window')
	# BITCoPTWindows = ['BITCoPT South Window']

print '_________________________________________________________________\n'

print "Modifying " + str(len(BITCoPTWindows)) + " window(s) for co-simulation with the BITCoPT system"

# Grab all the window objects in the file
windows = idf1.idfobjects['FENESTRATIONSURFACE:DETAILED']
# For each BITCoPT window being studied in this co-simulation
for bWindow in BITCoPTWindows:
	# For each window in the list of windows iteration though
	for window in windows:
		# If name of window matches the co-sim window then proceed 
		if window.Name == bWindow:

			print window

			# Add the SolarIncidentInside object to define 
			# solar heat gain in EnergyPlus during co-simulation
			idf1.newidfobject('SURFACEPROPERTY:SOLARINCIDENTINSIDE')
			solarInsideObj = idf1.idfobjects['SURFACEPROPERTY:SOLARINCIDENTINSIDE'][-1]
			solarInsideObj.Name = window.Name + ' Solar Incident Inside'
			solarInsideObj.Surface_Name = window.Name
			solarInsideObj.Construction_Name = window.Construction_Name
			solarInsideObj.Inside_Surface_Incident_Sun_Solar_Radiation_Schedule_Name = window.Name + ' SHG'

			print solarInsideObj

			# Add new object for the solar heat gain schedule 
			idf1.newidfobject('EXTERNALINTERFACE:FUNCTIONALMOCKUPUNITEXPORT:TO:SCHEDULE')
			toScheduleLast = idf1.idfobjects['EXTERNALINTERFACE:FUNCTIONALMOCKUPUNITEXPORT:TO:SCHEDULE'][-1]
			toScheduleLast.Schedule_Name = window.Name + ' SHG'
			toScheduleLast.Schedule_Type_Limits_Names = 'Any Number'
			toScheduleLast.FMU_Variable_Name = (window.Name + ' Solar Inside').replace(' ', '')
			toScheduleLast.Initial_Value = 0

			print toScheduleLast

			# Add new object for the cavity air temperature the respresents the 
			# new boundary condition in EnergyPlus for the BITCoPT cavity condition
			idf1.newidfobject('EXTERNALINTERFACE:FUNCTIONALMOCKUPUNITEXPORT:TO:ACTUATOR')
			toActuatorLast = idf1.idfobjects['EXTERNALINTERFACE:FUNCTIONALMOCKUPUNITEXPORT:TO:ACTUATOR'][-1]
			toActuatorLast.Name = (window.Name + ' Cavity Temp').replace(' ','')
			toActuatorLast.Actuated_Component_Unique_Name = window.Name
			toActuatorLast.Actuated_Component_Type = 'Surface'
			toActuatorLast.Actuated_Component_Control_Type = 'Outdoor Air Drybulb Temperature'
			toActuatorLast.FMU_Variable_Name = (window.Name + 'CavAirTemp').replace(' ', '')
			toActuatorLast.Initial_Value = 20

			print toActuatorLast

# Setup componets for daylighting input during co-simulation
lights = idf1.idfobjects['LIGHTS']
for direction in directions:
	for light in lights:
		if str(direction).lower() in (light.Name).lower():
			# Create new To:Schedule component that accepts the lighting fraction schedule value during co-simulation
			idf1.newidfobject('EXTERNALINTERFACE:FUNCTIONALMOCKUPUNITEXPORT:TO:SCHEDULE')
			toScheduleLast = idf1.idfobjects['EXTERNALINTERFACE:FUNCTIONALMOCKUPUNITEXPORT:TO:SCHEDULE'][-1]
			toScheduleLast.Schedule_Name = light.Name + ' Fractional Schedule'
			toScheduleLast.Schedule_Type_Limits_Names = 'Any Number'
			toScheduleLast.FMU_Variable_Name = (light.Name + 'Fraction').replace(' ', '')
			toScheduleLast.Initial_Value = 1
			# Set Lights component to use the new schedule name
			light.Schedule_Name = toScheduleLast.Schedule_Name
			light.Fraction_Replaceable = 1


# Setup components and schedules to connect BITCoPT hot water to 
# PlantComponent:TemperatureSource, Pump, and SetpointManager:Scheduled
bitcoptPlantComponents = []
for name in directions:
	bitcoptPlantComponents.append('BITCoPT ' + str(name) + ' Plant Component Temperature Source')

print '_________________________________________________________________\n'

print "Modifying "+str(len(bitcoptPlantComponents))+" plant loop(s) for co-simulation:"

setpointSchedules = idf1.idfobjects['SETPOINTMANAGER:SCHEDULED']
plantComponents = idf1.idfobjects['PLANTCOMPONENT:TEMPERATURESOURCE']
for direction in directions:
	for bitcopt in bitcoptPlantComponents:
		for component in plantComponents:
			if component.Name == bitcopt:
				# Add schedule input object to take BITCoPT fluid temperature
				idf1.newidfobject('EXTERNALINTERFACE:FUNCTIONALMOCKUPUNITEXPORT:TO:SCHEDULE')
				toScheduleLast = idf1.idfobjects['EXTERNALINTERFACE:FUNCTIONALMOCKUPUNITEXPORT:TO:SCHEDULE'][-1]
				toScheduleLast.Schedule_Name = component.Name + ' Schedule'
				toScheduleLast.Schedule_Type_Limits_Names = 'Any Number'
				toScheduleLast.FMU_Variable_Name = (str(direction) + 'FluidTemperature').replace(' ', '')
				toScheduleLast.Initial_Value = 0

				print toScheduleLast

				# Set PlantComponent:TemperatureSource scheduled temperature to the new input schedule
				component.Source_Temperature_Schedule_Name = toScheduleLast.Schedule_Name

				print component

				# Add actuator input object to control pump flowrate
				idf1.newidfobject('EXTERNALINTERFACE:FUNCTIONALMOCKUPUNITEXPORT:TO:ACTUATOR')
				toActuatorLast = idf1.idfobjects['EXTERNALINTERFACE:FUNCTIONALMOCKUPUNITEXPORT:TO:ACTUATOR'][-1]
				toActuatorLast.Name = ('BITCoPT '+str(direction)+' Var Spd Pump Mass Flow Actuator').replace(' ','')
				toActuatorLast.Actuated_Component_Unique_Name = 'BITCoPT '+str(direction)+' Var Spd Pump'
				toActuatorLast.Actuated_Component_Type = 'Pump'
				toActuatorLast.Actuated_Component_Control_Type = 'Pump Mass Flow Rate'
				toActuatorLast.FMU_Variable_Name = (str(direction)+'FluidFlow_kgps').replace(' ', '')
				toActuatorLast.Initial_Value = 0

				print toActuatorLast

				# Update setpoint manager schedule
				for schedule in setpointSchedules:
					if schedule.Name == 'BITCoPT '+str(direction)+' Scheduled HW Temp':
						schedule.Schedule_Name = toScheduleLast.Schedule_Name

						print schedule

print '_________________________________________________________________\n'

print "Setup generator and load center for BITCoPT electricity."

# Add generator components
turbine = idf1.idfobjects['GENERATOR:COMBUSTIONTURBINE']
if len(turbine) < 1:
	idf1.newidfobject('GENERATOR:COMBUSTIONTURBINE')
	turbine = idf1.idfobjects['GENERATOR:COMBUSTIONTURBINE'][-1]
	turbine.Name = 'BITCoPT Electricity Generator'
	turbine.Rated_Power_Output = 100000000.0 # Watts, this value is really just a maximum, as long as the number is high enough
	turbine.Minimum_Part_Load_Ratio = 0.0
	turbine.Maximum_Part_Load_Ratio = 1.0
	turbine.Optimum_Part_Load_Ratio = 0.5
	turbine.Part_Load_Based_Fuel_Input_Curve_Name = turbine.Name + ' Quadratic Curve'
	turbine.Temperature_Based_Fuel_Input_Curve_Name = turbine.Name + ' Quadratic Curve'
	turbine.Exhaust_Flow_Curve_Name = turbine.Name + ' Quadratic Curve'
	turbine.Part_Load_Based_Exhaust_Temperature_Curve_Name = turbine.Name + ' Quadratic Curve'
	turbine.Temperature_Based_Exhaust_Temperature_Curve_Name = turbine.Name + ' Quadratic Curve' 
	turbine.Heat_Recovery_Lube_Energy_Curve_Name = turbine.Name + ' Quadratic Curve'
	turbine.Coefficient_1_of_UFactor_Times_Area_Curve = 0.02
	turbine.Coefficient_2_of_UFactor_Times_Area_Curve = 0.9
	turbine.Maximum_Exhaust_Flow_per_Unit_of_Power_Output = 0.0
	turbine.Design_Minimum_Exhaust_Temperature = 150
	turbine.Design_Air_Inlet_Temperature = 25
	turbine.Fuel_Higher_Heating_Value = 43500
	turbine.Design_Heat_Recovery_Water_Flow_Rate = 0.0
	turbine.Fuel_Type = 'NaturalGas'

	idf1.newidfobject('GENERATOR:FUELSUPPLY')
	fuel = idf1.idfobjects['GENERATOR:FUELSUPPLY'][-1]
	fuel.Name = 'NaturalGas'
	fuel.Fuel_Temperature_Modeling_Mode = 'TemperatureFromAirNode'
	fuel.Fuel_Temperature_Reference_Node_Name = 'MicroCoGen1 air inlet node'
	fuel.Fuel_Temperature_Schedule_Name = ''                        
	fuel.Compressor_Power_Multiplier_Function_of_Fuel_Rate_Curve_Name = 'NullCubic'
	fuel.Compressor_Heat_Loss_Factor = 1.0000
	fuel.Fuel_Type = 'GaseousConstituents'
	fuel.Liquid_Generic_Fuel_Lower_Heating_Value = ''
	fuel.Liquid_Generic_Fuel_Higher_Heating_Value = ''
	fuel.Liquid_Generic_Fuel_Molecular_Weight = ''
	fuel.Liquid_Generic_Fuel_CO2_Emission_Factor = ''
	fuel.Number_of_Constituents_in_Gaseous_Constituent_Fuel_Supply = 8
	fuel.Constituent_1_Name = 'METHANE'
	fuel.Constituent_1_Molar_Fraction = 0.9490
	fuel.Constituent_2_Name = 'CarbonDioxide'
	fuel.Constituent_2_Molar_Fraction = 0.0070
	fuel.Constituent_3_Name = 'NITROGEN'
	fuel.Constituent_3_Molar_Fraction = 0.0160
	fuel.Constituent_4_Name = 'ETHANE'
	fuel.Constituent_4_Molar_Fraction = 0.0250
	fuel.Constituent_5_Name = 'PROPANE'
	fuel.Constituent_5_Molar_Fraction = 0.0020
	fuel.Constituent_6_Name = 'BUTANE'
	fuel.Constituent_6_Molar_Fraction = 0.0006
	fuel.Constituent_7_Name = 'PENTANE'
	fuel.Constituent_7_Molar_Fraction = 0.0002
	fuel.Constituent_8_Name = 'OXYGEN'
	fuel.Constituent_8_Molar_Fraction = 0.0002

	idf1.newidfobject('ELECTRICLOADCENTER:GENERATORS')
	generators = idf1.idfobjects['ELECTRICLOADCENTER:GENERATORS'][-1]
	generators.Name = turbine.Name + ' List'
	generators.Generator_1_Name = turbine.Name
	generators.Generator_1_Object_Type = 'Generator:CombustionTurbine'
	generators.Generator_1_Rated_Electric_Power_Output = turbine.Rated_Power_Output
	generators.Generator_1_Availability_Schedule_Name = ''
	generators.Generator_1_Rated_Thermal_to_Electrical_Power_Ratio = 1.0

	idf1.newidfobject('ELECTRICLOADCENTER:DISTRIBUTION')
	distribution = idf1.idfobjects['ELECTRICLOADCENTER:DISTRIBUTION'][-1]
	distribution.Name = turbine.Name + ' Distribution'
	distribution.Generator_List_Name = generators.Name
	distribution.Generator_Operation_Scheme_Type = 'TrackSchedule'
	distribution.Generator_Demand_Limit_Scheme_Purchased_Electric_Demand_Limit = 0
	distribution.Generator_Track_Schedule_Name_Scheme_Schedule_Name = turbine.Name + ' Schedule'
	distribution.Generator_Track_Meter_Scheme_Meter_Name = ''
	distribution.Electrical_Buss_Type = 'AlternatingCurrent'

	idf1.newidfobject('CURVE:QUADRATIC')
	quadratic = idf1.idfobjects['CURVE:QUADRATIC'][-1]
	quadratic.Name = turbine.Name + ' Quadratic Curve'
	quadratic.Coefficient1_Constant = 0.0
	quadratic.Coefficient2_x = 0.0
	quadratic.Coefficient3_x2 = 0.0
	quadratic.Minimum_Value_of_x = 0
	quadratic.Maximum_Value_of_x = 1

	idf1.newidfobject('CURVE:CUBIC')
	cubic = idf1.idfobjects['CURVE:CUBIC'][-1]
	cubic.Name = turbine.Name + ' Cubic Curve'
	cubic.Coefficient1_Constant = 0.0000
	cubic.Coefficient2_x = 0.0000
	cubic.Coefficient3_x2 = 0.0000
	cubic.Coefficient4_x3 = 0.0000
	cubic.Minimum_Value_of_x = 0.0000
	cubic.Maximum_Value_of_x = 0.0000

	print "Electricity generation from BITCoPT plugs into the schedule:"

	idf1.newidfobject('EXTERNALINTERFACE:FUNCTIONALMOCKUPUNITEXPORT:TO:SCHEDULE')
	toScheduleLast = idf1.idfobjects['EXTERNALINTERFACE:FUNCTIONALMOCKUPUNITEXPORT:TO:SCHEDULE'][-1]
	toScheduleLast.Schedule_Name = turbine.Name + ' Schedule'
	toScheduleLast.Schedule_Type_Limits_Names = 'Any Number'
	toScheduleLast.FMU_Variable_Name = 'BITCoPTElectricGeneration'
	toScheduleLast.Initial_Value = 0 

	print toScheduleLast





print '_________________________________________________________________\n'

print "Adding all desired output variables for studying co-simulation interactions, and setting reporting frequency to timestep"


# Add desired output variables to the file if they were not already a part of the file
outputVariables = [
	'Site Diffuse Solar Radiation Rate per Area',
	'Site Direct Solar Radiation Rate per Area',
	'Surface Outside Face Incident Solar Radiation Rate per Area',
	'Surface Outside Face Outdoor Air Drybulb Temperature',
	'Surface Outside Face Temperature',
	'Zone Thermostat Air Temperature',
	'Zone Total Internal Total Heating Energy',
	'Zone Windows Total Heat Gain Rate',
	'Zone Windows Total Heat Loss Rate',
	'Zone Windows Total Transmitted Solar Radiation Rate',
	'Zone Air System Sensible Cooling Rate',
	'Zone Air System Sensible Heating Rate',
	'Zone Air Temperature',
	'  Zone Mechanical Ventilation Mass Flow Rate',
	'Zone Infiltration Total Heat Loss Energy',
	'Zone Infiltration Total Heat Gain Energy',
	'Zone Lights Electric Energy',
	'Zone Lights Electric Power',
	'Zone Mean Air Temperature',
	'Zone Mean Radiant Temperature',
	'Zone Mean Air Humidity Ratio',
	'Zone Mechanical Ventilation Air Changes per Hour',
	'People Occupant Count',
	' Schedule Value',
	'Zone Mean Air Temperature',
	'Zone Outdoor Air Drybulb Temperature',
	'Schedule Value',
	'Surface Inside Face Absorbed Shortwave Radiation Rate',
	'Generator Produced Electric Power',
	'Generator Produced Electric Energy',
	'Generator Produced Thermal Rate',
	'Generator Produced Thermal Energy',
	'Electric Load Center Produced Electric Power',
	'Electric Load Center Produced Electric Energy',
	'Electric Load Center Supplied Electric Power',
	'Electric Load Center Requested Electric Power',
	'Plant Temperature Source Component Mass Flow Rate',
	'Plant Temperature Source Component Inlet Temperature',
	'Plant Temperature Source Component Outlet Temperature',
	'Plant Temperature Source Component Source Temperature',
	'Plant Temperature Source Component Heat Transfer Rate',
	'Plant Temperature Source Component Heat Transfer Energy',
	'Water Use Equipment Hot Water Mass Flow Rate',
	'Water Use Equipment Cold Water Mass Flow Rate',
	'Water Use Equipment Total Mass Flow Rate',
	'Water Use Equipment Hot Water Temperature',
	'Water Use Equipment Cold Water Temperature',
	'Water Use Equipment Target Water Temperature',
	'Water Use Equipment Mixed Water Temperature',
	'Water Use Equipment Total Volume',
	'Fluid Heat Exchanger Heat Transfer Rate',
	'Fluid Heat Exchanger Heat Transfer Energy',
	'Fluid Heat Exchanger Loop Supply Side Mass Flow Rate',
	'Fluid Heat Exchanger Loop Supply Side Inlet Temperature',
	'Fluid Heat Exchanger Loop Supply Side Outlet Temperature',
	'Fluid Heat Exchanger Loop Demand Side Mass Flow Rate',
	'Fluid Heat Exchanger Loop Demand Side Inlet Temperature',
	'Fluid Heat Exchanger Loop Demand Side Outlet Temperature',
	'Fluid Heat Exchanger Operation Status',
	'Fluid Heat Exchanger Effectiveness',
	'Baseboard Total Heating Rate',
	'Baseboard Hot Water Mass Flow Rate',
	'Baseboard Water Inlet Temperature',
	'Baseboard Water Outlet Temperature'
	]

# Get all output variables already in the file
allOutputs = idf1.idfobjects['OUTPUT:VARIABLE']

# Check if there are any output variables in the IDF
# If no output variables, just add what we want
if len(allOutputs) < 1:
	for variable in outputVariables:
		idf1.newidfobject('OUTPUT:VARIABLE')
		lastOutput = idf1.idfobjects['OUTPUT:VARIABLE'][-1]
		lastOutput.Variable_Name = str(variable)
		lastOutput.Key_Value = '*'
		lastOutput.Reporting_Frequency = 'Timestep'
# else check the we don't have duplicate output variables
# and add the desired output variables
else:
	for variable in outputVariables:
		for output in allOutputs:
			# Check if there is already an output variable for what we want
			if output.Variable_Name == variable:
				# Update output variable for Timestep reporting frequency
				output.Reporting_Frequency = 'Timestep'
			# Add all the remaining desired output variables
			else:
				idf1.newidfobject('OUTPUT:VARIABLE')
				lastOutput = idf1.idfobjects['OUTPUT:VARIABLE'][-1]
				lastOutput.Variable_Name = str(variable)
				lastOutput.Key_Value = '*'
				lastOutput.Reporting_Frequency = 'Timestep'

# print 'All output variables:'
# print allOutputs

SQLite = idf1.idfobjects['OUTPUT:SQLITE'][-1]
SQLite.Option_Type = ''

print '_________________________________________________________________\n'

newFile = fname.replace('.idf','_CoSim.idf')

print 'Saving IDF file at: '+str(newFile)+'\n'

idf1.saveas(newFile)

print 'The IDF file is now ready to be packaged as an FMU through FMU Export using the EnergyPlusToFMU script developed by LBNL at: http://simulationresearch.lbl.gov/fmu/EnergyPlus/export/index.html'