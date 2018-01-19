#!/usr/bin/env python
import math,struct,os,re
import numpy as np
from shell_data import *
from trimesh_data_structure import *
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
    ############################################################################
    ## @param quat quaternion that is vector-4 [w, x, y, z], w is the real part
    ############################################################################
    @staticmethod
    def Quaternion_To_Axis_Rotation_Degree(quat):
        rRad = 2.0 * np.arccos(quat[0])
        if (rRad < 1E-12):
            rAxis = np.array([0., 1., 0.])
        else:
            rAxis = np.array(quat[1:])
            rAxis = rAxis / np.linalg.norm(rAxis)
        rDeg = 180./np.pi * rRad
        return rAxis, rDeg
    @staticmethod
    def Quaternion_Slerp(quat1, quat2, alpha):
        '''Referenced: https://stackoverflow.com/questions/2879441/how-to-interpolate-rotations and
         https://en.wikipedia.org/wiki/Slerp '''
        p0, p1 = np.array(quat1), np.array(quat2)
        dot = np.dot(p0 / np.linalg.norm(p0), p1 / np.linalg.norm(p1))
        if dot < 0.0: # opposite handedness so must reverse so slerp takes shorter path
            p1 = -p1
            dot = -dot
        DOT_THRESH = 0.9995
        if dot > DOT_THRESH: # too close just, just Lerp
            new_quat = p0*(1.0 - alpha) + p1*alpha
            new_quat = new_quat / np.linalg.norm(new_quat)
            return new_quat.tolist()
        omega = np.arccos(dot)
        so = np.sin(omega)
        new_quat = (np.sin((1.0 - alpha) * omega) / so) * p0 + (np.sin(alpha * omega) / so) * p1
        return new_quat.tolist()

################################################################################
################################################################################
class Translation:
    def __init__(self, translation):
        assert(len(translation)==3)
        self.values = translation[:]
    def Format_String(self):
        s = """            <translate x="%f" y="%f" z="%f"/>""" %(self.values[0], self.values[1], self.values[2])
        return s
    @staticmethod
    def Lerp(trans1, trans2, alpha):
        ''' Linear interpolation between 2 translations. '''
        interp_trans = np.array(trans1)*(1.0 - alpha) + np.array(trans2)*alpha
        return interp_trans.tolist()

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
class Object(object):
    def __init__(self):
        self.restCOM = None
        self.ID = -1
        self.material_format_string = None
        print 'run Object class constructor'
    def Read_Rest_COM(self, filename):
        assert(os.path.isfile(filename))
        ifs = open(filename, 'rb')
        all_bytes = ifs.read()
        self.restCOM = struct.unpack('ddd', all_bytes)
    def Set_Material_String(self, s_mat):
        self.material_format_string = s_mat

################################################################################
################################################################################
class Rigid_Wavefront_Obj(Object):
    def __init__(self, filename, use_COM=True):
        super(Rigid_Wavefront_Obj, self).__init__()
        self.objname = None
        self.filename = filename
        if use_COM:
            base = os.path.basename(filename)
            self.objname = base.split('.')[0]
            prefix = os.path.dirname(filename)
            print prefix, self.objname
            rest_com_file = '%s/%s_centerOfMass.3vector' %(prefix, self.objname)
            self.Read_Rest_COM(rest_com_file)
        self.frames = []
    def Format_String(self, frame, material_string=''):
        if (self.material_format_string is not None):
            material_string = self.material_format_string
        if (frame != -1):
            frame_string = self.frames[frame].Format_String()
        else:
            frame_string = ''
        s = """
    <shape type="obj">
        <string name="filename" value="%s"/> %s %s
    </shape> """ %(self.filename, frame_string, material_string)
        return s

################################################################################
################################################################################
class Ply(Object):
    def __init__(self, filename):
        self.filename = filename
        self.frames = []
    def Format_String(self, frame, material_string=''):
        if (self.material_format_string is not None):
            material_string = self.material_format_string
        s = """
    <shape type="ply">
        <string name="filename" value="%s"/> %s %s
        <bsdf type="diffuse">
         <texture type="vertexcolors" name="reflectance"/>
        </bsdf>
    </shape> """ %(self.filename, self.frames[frame].Format_String(), material_string)
        return s

################################################################################
################################################################################
class Shell_Obj(Object):
    def __init__(self, filename):
        super(Shell_Obj, self).__init__()
        self.filename = filename
        base = os.path.basename(filename)
        self.objname = base.split('.')[0]
        prefix = os.path.dirname(filename)
        self.frames = dict()
        self.v_map = None
        self.match_re = None
        self.render2sim_bary_map = []
        self.sim_obj = None
        self.map_read = False
    def Set_MatchRegex(self, re):
        self.match_re = re
    def Read_Data(self, directory, render2sim_file):
        if not self.map_read:
            Read_Shell_Vertex_Map(directory, self)
            Read_Shell_Render2Sim_Bary_Map(self, render2sim_file)
        Read_Shell_Disp(directory, self)
        self.map_read = True
    def Format_String(self, frame, material_string=''):
        if (self.material_format_string is not None):
            material_string = self.material_format_string

        # Figure out the deformation first
        deformed_obj = TriangleMesh()
        deformed_obj.Read(self.filename)
        sim_obj = TriangleMesh()
        sim_obj.Read('../assets/cymbal/cymbal_remeshed_18inch.obj')

        # update underlying sim obj vertices
        N_udofs = len(frame.disp)
        N_uvtxs = N_udofs/3
        # uses this to hold displacements
        sim_obj.v = np.zeros(sim_obj.v.shape)
        for ii in range(N_uvtxs):
            original_v = self.v_map[ii,1]
            sim_obj.v[original_v,0] = frame.disp[ii*3+0]
            sim_obj.v[original_v,1] = frame.disp[ii*3+1]
            sim_obj.v[original_v,2] = frame.disp[ii*3+2]

        ii = 0
        for bary in self.render2sim_bary_map:
            vidx0_o = sim_obj.f[bary[0],0]
            vidx1_o = sim_obj.f[bary[0],1]
            vidx2_o = sim_obj.f[bary[0],2]
            dis = bary[1][0]*sim_obj.v[vidx0_o, :] \
                + bary[1][1]*sim_obj.v[vidx1_o, :] \
                + bary[1][2]*sim_obj.v[vidx2_o, :]
            deformed_obj.v[ii,:] += dis
            ii += 1

        deformed_obj_filename = 'tmp.obj'
        deformed_obj.Write(deformed_obj_filename)

        s = """
    <shape type="obj">
        <string name="filename" value="%s"/> %s
    </shape> """ %(deformed_obj_filename, material_string)
        return s

if __name__ == "__main__":
    # test
    obj = Rigid_Wavefront_Obj('/home/jui-hsien/data/models/B/B.obj')
    obj.Add_Rotation([1.,0.,0.], 180.)
    obj.Add_Translation([2.,3.,4.])
    print obj.Format_String()
