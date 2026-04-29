import adsk.core
import adsk.fusion
import os
import csv
from ...lib import fusion360utils as futil
from ... import config
app = adsk.core.Application.get()
ui = app.userInterface


CMD_ID = f'{config.COMPANY_NAME}_{config.ADDIN_NAME}_PropellerBuilder'
CMD_NAME = 'Propeller Builder'
CMD_Description = 'Loads control or fit points into sketches. Lofts into propeller. Option for STEP export.'
IS_PROMOTED = True

# TODO *** Define the location where the command button will be created. ***
# This is done by specifying the workspace, the tab, and the panel, and the 
# command it will be inserted beside. Not providing the command to position it
# will insert it at the end.
WORKSPACE_ID = 'FusionSolidEnvironment'
PANEL_ID = 'SolidScriptsAddinsPanel'
COMMAND_BESIDE_ID = 'ScriptsManagerCommand'

# Resource location for command icons, here we assume a sub folder in this directory named "resources".
ICON_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'resources', '')

# Local list of event handlers used to maintain a reference so
# they are not released and garbage collected.
local_handlers = []


# Executed when add-in is run.
def start():
    # Create a command Definition.
    cmd_def = ui.commandDefinitions.addButtonDefinition(CMD_ID, CMD_NAME, CMD_Description, ICON_FOLDER)

    # Define an event handler for the command created event. It will be called when the button is clicked.
    futil.add_handler(cmd_def.commandCreated, command_created)

    # ******** Add a button into the UI so the user can run the command. ********
    # Get the target workspace the button will be created in.
    workspace = ui.workspaces.itemById(WORKSPACE_ID)

    # Get the panel the button will be created in.
    panel = workspace.toolbarPanels.itemById(PANEL_ID)

    # Create the button command control in the UI after the specified existing command.
    control = panel.controls.addCommand(cmd_def, COMMAND_BESIDE_ID, False)

    # Specify if the command is promoted to the main toolbar. 
    control.isPromoted = IS_PROMOTED


# Executed when add-in is stopped.
def stop():
    # Get the various UI elements for this command
    workspace = ui.workspaces.itemById(WORKSPACE_ID)
    panel = workspace.toolbarPanels.itemById(PANEL_ID)
    command_control = panel.controls.itemById(CMD_ID)
    command_definition = ui.commandDefinitions.itemById(CMD_ID)

    # Delete the button command control
    if command_control:
        command_control.deleteMe()

    # Delete the command definition
    if command_definition:
        command_definition.deleteMe()


# Function that is called when a user clicks the corresponding button in the UI.
# This defines the contents of the command dialog and connects to the command related events.
def command_created(args: adsk.core.CommandCreatedEventArgs):
    # General logging for debug.
    futil.log(f'{CMD_NAME} Command Created Event')

    # https://help.autodesk.com/view/fusion360/ENU/?contextId=CommandInputs
    inputs = args.command.commandInputs

    # Create the button for file in selection.
    inputs.addBoolValueInput('fileInID', 'Select Input File', False)

    # Create the file in selection box.
    inputs.addTextBoxCommandInput('fileInBoxID', 'Input File ID', "", 1, False)

    # Create the radio button toggle for control/fit.
    pointTypeID = inputs.addRadioButtonGroupCommandInput('pointTypeID', 'Point Type')
    radioButtonItems = pointTypeID.listItems
    radioButtonItems.add("Control Point Spline", True)
    radioButtonItems.add("Fit Point Spline", False)

    # Create the boolean toggle for step file output.
    inputs.addBoolValueInput('stepID', 'STEP Export?', True, "", True)
    
    # Create the boolean toggle for fixed splines.
    inputs.addBoolValueInput('fixedID', 'Fix Splines?', True, "", True)

    # Create the button for folder out selection.
    inputs.addBoolValueInput('folderOutID', 'Select Output Folder', False)

    # Create the folder out selection box.
    inputs.addTextBoxCommandInput('folderOutBoxID', 'Output Folder ID', "", 1, False)     
    
    # Create the file out selection box.
    inputs.addTextBoxCommandInput('fileOutID', 'Output File ID', "", 1, False)    


    # TODO Connect to the events that are needed by this command.
    futil.add_handler(args.command.execute, command_execute, local_handlers=local_handlers)
    futil.add_handler(args.command.inputChanged, command_input_changed, local_handlers=local_handlers)
    futil.add_handler(args.command.executePreview, command_preview, local_handlers=local_handlers)
    futil.add_handler(args.command.validateInputs, command_validate_input, local_handlers=local_handlers)
    futil.add_handler(args.command.destroy, command_destroy, local_handlers=local_handlers)


# This event handler is called when the user clicks the OK button in the command dialog or 
# is immediately called after the created event not command inputs were created for the dialog.
    
#Jules - Does the thing
def command_execute(args: adsk.core.CommandEventArgs):
    # General logging for debug.
    futil.log(f'{CMD_NAME} Command Execute Event')

    # Get a reference to your command's inputs.
    inputs = args.command.commandInputs
    fileInID: adsk.core.BoolValueCommandInput = inputs.itemById('fileInID')
    fileInBoxID: adsk.core.TextBoxCommandInput = inputs.itemById('fileInBoxID')
    pointTypeID: adsk.core.RadioButtonGroupCommandInput = inputs.itemById('pointTypeID')
    fixedID: adsk.core.BoolValueCommandInput = inputs.itemById('fixedID')
    stepID: adsk.core.BoolValueCommandInput = inputs.itemById('stepID')
    folderOutID: adsk.core.BoolValueCommandInput = inputs.itemById('folderOutID')
    folderOutBoxID: adsk.core.TextBoxCommandInput = inputs.itemById('folderOutBoxID')
    fileOutID: adsk.core.TextBoxCommandInput = inputs.itemById('fileOutID')
    
    rootComp = app.activeProduct.rootComponent
    sketches = rootComp.sketches 
    guideCurve = False 
    density = 0
    bladeCount = 0

    #Jules - reads in points, does different things for each
    with open(fileInBoxID.text, "r") as f:
        reader = csv.reader(f, delimiter = "\t")
        for row in reader:
            if (row[0][0:3] == "Air"): #makes a new sketch
                sketch = sketches.add(rootComp.xYConstructionPlane)
                sketch.name = row[0]
                guideCurve = False
                density += 1
            elif (row[0][0:3] == "Gui"): #makes a new sketch
                sketch = sketches.add(rootComp.xYConstructionPlane)
                sketch.name = row[0]
                guideCurve = True
            elif (row[0][0:3] == "Hub"): #makes a new sketch with a cylinder
                sketch = sketches.add(rootComp.xYConstructionPlane)
                substrings = str.split(row[0])
                sketch.name = substrings[0]
                sketch.sketchCurves.sketchCircles.addByCenterRadius(adsk.core.Point3D.create(0, 0, 0), float(substrings[1]))
            elif (row[0][0:3] == "Bla"): #sets the blade number
                substrings = str.split(row[0])
                bladeCount = substrings[1]
            elif (row[0] == "START"): #initializer
                controlPoints = []
                fitPoints = adsk.core.ObjectCollection.create()
            elif (row[0] == "END"): #publishes points to sketches
                if (guideCurve):
                    sketch.sketchCurves.sketchFittedSplines.add(fitPoints)
                else:
                    if (pointTypeID.selectedItem.index == 0):
                        deg = 0
                        if (len(controlPoints) >= 6):
                            deg = 5
                        else:
                            deg = len(controlPoints)
                        sketch.sketchCurves.sketchControlPointSplines.add(controlPoints, deg)
                    else:
                        sketch.sketchCurves.sketchFittedSplines.add(fitPoints)
            else: #adds points to arrays
                controlPoints.append(adsk.core.Point3D.create(float(row[0]), float(row[1]), float(row[2])))
                fitPoints.add(adsk.core.Point3D.create(float(row[0]), float(row[1]), float(row[2])))
    
    #Jules - fixes sketches if necessary
    if fixedID.value:
        for s in sketches:
            for c in s.sketchCurves:
                c.isFixed = True
    
    #Jules - surface lofting and hub extrusion   
    loftFeats = rootComp.features.loftFeatures
    loftInput = loftFeats.createInput(adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
    loftSectionsObj = loftInput.loftSections
    loftRailsObj = loftInput.centerLineOrRails
    
    for s in sketches: 
        if (s.name[0:3] == "Air"):
            airfoil = adsk.core.ObjectCollection.create()
            for cps in s.sketchCurves.sketchControlPointSplines:
                airfoil.add(cps)
            for sfs in s.sketchCurves.sketchFittedSplines:
                airfoil.add(sfs)
            loftSectionsObj.add(adsk.fusion.Path.create(airfoil, adsk.fusion.ChainedCurveOptions.tangentChainedCurves))
        elif (s.name[0:3] == "Gui"):
            rail = s.sketchCurves.item(0)
            loftRailsObj.addRail(rail)
        elif (s.name[0:3] == "Hub"):
            extInput = rootComp.features.extrudeFeatures.createInput(s.profiles.item(0), adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
            extInput.setSymmetricExtent(adsk.core.ValueInput.createByReal(3), True) #hub length parameter
            rootComp.features.extrudeFeatures.add(extInput)

    loftInput.isSolid = False
    loftInput.isClosed = False
    loftInput.isTangentEdgesMerged = True
    loftFeats.add(loftInput)

    #Jules - end patching
    patchFeats = rootComp.features.patchFeatures
    airfoil1 = adsk.core.ObjectCollection.create()
    for cps in sketches.itemByName("Airfoil1").sketchCurves.sketchControlPointSplines:
        airfoil1.add(cps)
    for sfs in sketches.itemByName("Airfoil1").sketchCurves.sketchFittedSplines:
        airfoil1.add(sfs)
    patchPath1 = adsk.fusion.Path.create(airfoil1, adsk.fusion.ChainedCurveOptions.tangentChainedCurves)
    patchInput1 = patchFeats.createInput(patchPath1, adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
    patchFeats.add(patchInput1)

    airfoil2 = adsk.core.ObjectCollection.create()
    for cps in sketches.item(density - 1).sketchCurves.sketchControlPointSplines:
        airfoil2.add(cps)
    for sfs in sketches.item(density - 1).sketchCurves.sketchFittedSplines:
        airfoil2.add(sfs)
    patchPath2 = adsk.fusion.Path.create(airfoil2, adsk.fusion.ChainedCurveOptions.tangentChainedCurves)
    patchInput2 = patchFeats.createInput(patchPath2, adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
    patchFeats.add(patchInput2)
    
    #Jules - stitches patches to loft
    stitchBodies = adsk.core.ObjectCollection.create()
    bodies = rootComp.bRepBodies
    stitchBodies.add(bodies.item(1))
    stitchBodies.add(bodies.item(2))
    stitchBodies.add(bodies.item(3))
    stitchFeats = rootComp.features.stitchFeatures
    tolerance = adsk.core.ValueInput.createByReal(.00001)
    stitchInput = stitchFeats.createInput(stitchBodies, tolerance)
    stitchFeats.add(stitchInput)
    
    #Jules - patterns the blades around the z-axis
    blade = adsk.core.ObjectCollection.create()
    bodies = rootComp.bRepBodies
    blade.add(bodies.item(1))
    circPattFeat = rootComp.features.circularPatternFeatures
    zAxis = rootComp.zConstructionAxis
    circInput = circPattFeat.createInput(blade, zAxis)
    circInput.quantity = adsk.core.ValueInput.createByReal(float(bladeCount))
    circInput.totalAngle = adsk.core.ValueInput.createByString('360 deg')
    circInput.isSymmetric = False
    circPattFeat.add(circInput)
    
    #Jules - combines everything
    combineBodies = adsk.core.ObjectCollection.create()
    bodies = rootComp.bRepBodies
    for b in range(1, bodies.count):
        combineBodies.add(bodies.item(b))    
    combineFeats = rootComp.features.combineFeatures
    combineInput = combineFeats.createInput(bodies.item(0), combineBodies)
    combineFeats.add(combineInput)

    #Jules - STEP exporting
    if stepID.value:
        exportMgr = app.activeProduct.exportManager
        stepOptions = exportMgr.createSTEPExportOptions(folderOutBoxID.text + "/" + fileOutID.text)
        exportMgr.execute(stepOptions)

# This event handler is called when the command needs to compute a new preview in the graphics window.
def command_preview(args: adsk.core.CommandEventArgs):
    # General logging for debug.
    futil.log(f'{CMD_NAME} Command Preview Event')
    inputs = args.command.commandInputs


# This event handler is called when the user changes anything in the command dialog
# allowing you to modify values of other inputs based on that change.

#Jules - defines what happens when either button is clicked
def command_input_changed(args: adsk.core.InputChangedEventArgs):
    changed_input = args.input
    inputs = args.inputs

    if changed_input.id == 'fileInID':
        file_dialog = ui.createFileDialog()
        file_dialog.title = "Select a file"
        file_dialog.filter = 'Text files (*.txt)'
        file_dialog.filterIndex = 0
        file_dialog.isMultiSelectEnabled = False
        file_dialog.showOpen()
        inputs.itemById('fileInBoxID').text = str(file_dialog.filename)

    if changed_input.id == 'folderOutID':
        folder_dialog = ui.createFolderDialog()
        folder_dialog.title = "Select a folder"
        folder_dialog.showDialog()
        inputs.itemById('folderOutBoxID').text = str(folder_dialog.folder)

    # General logging for debug.
    futil.log(f'{CMD_NAME} Input Changed Event fired from a change to {changed_input.id}')


# This event handler is called when the user interacts with any of the inputs in the dialog
# which allows you to verify that all of the inputs are valid and enables the OK button.
def command_validate_input(args: adsk.core.ValidateInputsEventArgs):
    # General logging for debug.
    futil.log(f'{CMD_NAME} Validate Input Event')

    inputs = args.inputs
    
    # Verify the validity of the input values. This controls if the OK button is enabled or not.
    valueInput = inputs.itemById('value_input')
    if valueInput.value >= 0:
        args.areInputsValid = True
    else:
        args.areInputsValid = False
        

# This event handler is called when the command terminates.
def command_destroy(args: adsk.core.CommandEventArgs):
    # General logging for debug.
    futil.log(f'{CMD_NAME} Command Destroy Event')

    global local_handlers
    local_handlers = []
