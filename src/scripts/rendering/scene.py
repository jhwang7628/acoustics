#!/usr/bin/env python 
import math
################################################################################
################################################################################
class Rotation:
    def __init__(self):
        self.values = [0.0, 0.0, 0.0, 0.0] 
        self.type = None # either 'quaternion' or 'axis_angle_deg'
    def Set_Quaternion(self, values):
        assert(len(values)==4)
        self.values = values[:]
        self.type = 'quaternion'
    def Set_Axis_Angle_Deg(self, rotation_axis, rotation_angle_deg): 
        self.values = [rotation_axis[0], rotation_axis[1], rotation_axis[2], rotation_angle_deg]
        self.type = 'axis_angle_deg'
    def Is_Set(self):
        return self.type is not None
    def Format_String(self):
        assert(self.Is_Set())
        if self.type == 'axis_angle_deg': 
            # NOTE: it might be nan here for anxis, angle, manually check frames and decide what to do
            if math.isnan(self.values[3]): 
                self.values[3] = 0.0
            s = """            <rotate x="%f" y="%f" z="%f" angle="%f"/>""" %(self.values[0], self.values[1], self.values[2], self.values[3])
        else: 
            assert(False)
        return s

################################################################################
################################################################################
class Translation:
    def __init__(self, translation):
        assert(len(translation)==3)
        self.values = translation[:]
    def Format_String(self):
        s = """            <translate x="%f" y="%f" z="%f"/>""" %(self.values[0], self.values[1], self.values[2])
        return s

################################################################################
################################################################################
class Rigid_Frame: 
    def __init__(self, ID):
        self.ID = ID
        self.transforms = [] # first one in the list will be added to config first, and thus first applied
    def Add_Rotation(self, rotation_axis, rotation_angle_deg): 
        transform = Rotation()
        transform.Set_Axis_Angle_Deg(rotation_axis, rotation_angle_deg)
        self.transforms.append(transform)
    def Add_Translation(self, translation): 
        transform = Translation(translation)
        self.transforms.append(transform)
    def Format_String(self):
        string_transforms = '\n'
        if len(self.transforms) > 0:
            string_transforms += '        <transform name="toWorld">\n'
            for t in self.transforms:
                string_transforms += t.Format_String() + '\n'
            string_transforms += '        </transform>'
        else:
            string_transforms = ''
        return string_transforms
            
################################################################################
################################################################################
class Object: 
    def __init__(self): 
        self.ID = -1

################################################################################
################################################################################
class Rigid_Wavefront_Obj(Object): 
    def __init__(self, filename): 
        self.filename = filename
        self.frames = []
    def Format_String(self, frame, material_string=''):
        s = """
    <shape type="obj">
        <string name="filename" value="%s"/> %s %s
    </shape> """ %(self.filename, self.frames[frame].Format_String(), material_string)
        return s

################################################################################
################################################################################
class Ply(Object): 
    def __init__(self, filename): 
        self.filename = filename
        self.frames = []
    def Format_String(self, frame, material_string=''):
        s = """
    <shape type="ply">
        <string name="filename" value="%s"/> %s %s
        <bsdf type="diffuse">
         <texture type="vertexcolors" name="reflectance"/>
        </bsdf>
    </shape> """ %(self.filename, self.frames[frame].Format_String(), material_string)
        return s
       
if __name__ == "__main__":
    # test
    obj = Rigid_Wavefront_Obj('/home/jui-hsien/data/models/B/B.obj')
    obj.Add_Rotation([1.,0.,0.], 180.)
    obj.Add_Translation([2.,3.,4.])
    print obj.Format_String()
