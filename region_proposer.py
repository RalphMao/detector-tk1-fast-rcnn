import subprocess
import numpy as np
'''
#import sys
#sys.path.insert(0,'/home/maohz12/Downloads')
#import selective_search_ijcv_with_python 
#import os
#from PIL import Image
#script_dirname = os.path.abspath(os.path.dirname(__file__)) + '/'
selective_search = selective_search_ijcv_with_python.get_windows

class proposer(object):
    def __init__(self,max_num_proposals):
        self.max_num_proposals = max_num_proposals

    def get_proposals(self, image_fname):
        if image_fname[0] != '/':
            image_f = script_dirname + image_fname
        else:
            image_f = image_fname
        windows = selective_search([image_f])[0]
        bboxs =  map(lambda x:(x[1],x[0],x[3],x[2]), windows)
        if len(bboxs) < 20:
            image = Image.open(image_f)
            width, height = image.size
            slice_num = 4
            for i in range(slice_num):
                for j in range(slice_num):
                    bboxs.append((width*i/slice_num, height*j/slice_num, width * (i+1) / slice_num, height * (j+1)/slice_num))
        return bboxs
        
'''
'''
class proposer(object):
    def __init__(self,max_num_proposals):
        self.max_num_proposals = max_num_proposals

    def get_proposals(self, image_fname):
        image = Image.open(image_fname)
        width, height = image.size
        slice_num = 6
        bboxs = []
        for i in range(slice_num):
            for j in range(slice_num):
                bboxs.append((width*i/slice_num, height*j/slice_num, width * (i+1) / slice_num, height * (j+1)/slice_num))
        return bboxs
''' 
class proposer(object):
    def __init__(self,max_num_proposals):
        self.max_num_proposals = max_num_proposals

    def get_proposals(self, image_fname, dtype = 'float32'):
        string = subprocess.Popen(["./myedgeboxes", image_fname, str(self.max_num_proposals)],stdout = subprocess.PIPE).communicate()[0]
        bboxes = np.array(map(lambda x: int(x), string.split()), dtype = dtype)
        bboxes = np.reshape(bboxes, (len(bboxes) / 4, 4))
        # bboxes = bboxes[:,np.array([1,0,3,2])]
        return bboxes

if __name__ == "__main__":
    import os
    import PIL.Image as Image
    pp = proposer(200)
    print pp.get_proposals('edgeboxes-c/peppers.png')
    image = np.array(Image.open('edgeboxes-c/peppers.png'))
    print image.shape
