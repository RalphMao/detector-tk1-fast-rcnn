import subprocess
import numpy as np

class proposer(object):
    def __init__(self,max_num_proposals):
        self.max_num_proposals = max_num_proposals

    def get_proposals(self, image_fname, dtype = 'float32'):
        string = subprocess.Popen(["./myedgeboxes", image_fname, 'models/forest/modelBsds.dat'],stdout = subprocess.PIPE).communicate()[0]
        if string.strip() == '':
            return []
        bboxes = np.array(map(lambda x: int(x), string.split()), dtype = dtype)
        bboxes = np.reshape(bboxes, (len(bboxes) / 4, 4))
        # bboxes = bboxes[:,np.array([1,0,3,2])]
        return bboxes

if __name__ == "__main__":
    import os
    import sys
    pp = proposer(200)
    print pp.get_proposals(sys.argv[1])
