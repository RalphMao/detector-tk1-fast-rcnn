import subprocess
import time
import numpy as np
import PIL.Image as Image

class proposer(object):
    def __init__(self,max_num_proposals):
        self.max_num_proposals = max_num_proposals

    def get_proposals(self, image_fname, dtype = 'float32'):
        print image_fname
        string = subprocess.Popen(["./myedgeboxes", image_fname, 'models/forest/modelBsds.dat'],stdout = subprocess.PIPE).communicate()[0]
        if string.strip() == '':
            return []
        bboxes = np.array(map(lambda x: int(x), string.split()), dtype = dtype)
        bboxes = np.reshape(bboxes, (len(bboxes) / 4, 4))
        # bboxes = bboxes[:,np.array([1,0,3,2])]
        return bboxes[:self.max_num_proposals]

if __name__ == "__main__":
    import os
    import sys
    pp = proposer(200)
    Size = [300,400, 600, 800]
    for i_size in Size:
        image = Image.open(sys.argv[1])
        image_t = image.resize((i_size, i_size))
        image_t.save('images/image%d.jpg'%i_size)
        start_time = time.time()
        for x in range(10):
            pp.get_proposals('images/image%d.jpg'%i_size)
        print "Num %d, Time %f"%(i_size, time.time()-start_time) 

