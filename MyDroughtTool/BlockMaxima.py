import arcpy as ARCPY
import arcpyWithR as RARC
import numpy as NUM
import ErrorUtils as ERROR
import subprocess as SUB
import os as OS
import sys as SYS
import locale as LOCALE
LOCALE.setlocale(LOCALE.LC_ALL, '')

def BlockMaxima():
    #### Get User Provided Inputs ####
    inputWorkSpace = ARCPY.GetParameterAsText(0)    
    outputWorkSpace = ARCPY.GetParameterAsText(1) 
    rperiods = ARCPY.GetParameterAsText(2)
    rperiods = [str(i) for i in rperiods.split(";")]
    rperiods = ";".join(rperiods)
    useMinima = ARCPY.GetParameterAsText(3)
    useConfInt = ARCPY.GetParameterAsText(4) 
    useCoeff = ARCPY.GetParameterAsText(5)
	
    if useMinima == 'true':
        useMinima = "1"
    else:
        useMinima = "0"

    if useConfInt == 'true':
        useConfInt = "1"
    else:
        useConfInt = "0"
	
    if useCoeff == 'true':
        useCoeff = "1"
    else:
        useCoeff = "0"
						
    #### Create R Command ####
    pyScript = SYS.argv[0]
    toolDir = OS.path.dirname(pyScript)
    rScript = OS.path.join(toolDir, "BlockMaxima.r")
    ARCPY.SetProgressor("default", "Executing R Script...")

    args = ["R", "--slave", "--vanilla", "--args", rScript, inputWorkSpace, 
           outputWorkSpace, rperiods, useMinima, useConfInt, useCoeff]
		
    #### Uncomment Next Two Lines to Print/Create Command Line Args ####
    #md = RARC.createRCommand(args, rScript)
    #ARCPY.AddWarning(cmd)

    #### Execute Command ####
    scriptSource = open(rScript, 'rb')
    rCommand = SUB.Popen(args, 
                         stdin = scriptSource,
                         stdout = SUB.PIPE, 
                         stderr = SUB.PIPE,
                         shell=True)

    #### Print Result ####
    resString, errString = rCommand.communicate()

    #### Push Output to Message Window ####
    if errString and "Calculations Complete..." not in resString:
        ARCPY.AddError(errString)
    else:
        resOutString = RARC.printRMessages(resString)
        ARCPY.AddMessage(resOutString)            
	
if __name__ == '__main__':
    test = BlockMaxima() 