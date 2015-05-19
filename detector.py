'''
Author: Ralph Mao
May, 2015
'''
import edgebox
import nms
import cnn
import httpapi
from PIL import Image
import numpy as np

class rcnn_detector(object):
    def __init__(self, prototxt_name, model_name, batch_size = 1, max_num_proposals = 200, output_dim = 201, iou_thres = 0.3):
        self.Net = cnn.Net(prototxt_name, model_name, batch_size)
        self.edgebox = edgebox.detector(max_num_proposals)
        self.nms = nms.reducer(iou_thres)

    def detect(self, image_in):
        
        bboxs = self.edgebox.get_proposals(image_in)
        if type(image_in) == type(''): # if the path is given
            image = np.array(Image.load(image_in),dtype='float32')/255
        else: # if image data are given
            image = np.array(image_in,dtype='float32')/255
        regions = self.get_regions(image,bboxs)
        probs = np.zeros(regions.shape[0], self.output_dim, dtype = 'float32')
        batch_size = self.Net.max_batch
        for start in range(0, regions.shape[0], batch_size):
            np.copyto(probs[start:start + batch_size], self.Net.forward(regions[start:start + batch_size])
        keep = self.nms.nmsreduce(bboxs, probs[:,:-1]) # Other parameters?
        return keep

class http_rcnn_detector(rcnn_detector):
    def __init__(self, prototxt_name, model_name, batch_size = 1, max_num_proposals = 200, output_dim = 201, iou_thres = 0.3, usrname, passwd):
        rcnn_detector.__init__(prototxt_name, model_name, batch_size, max_num_proposals, output_dim, iou_thres)
        self.carrier = httpapi.carrier(usrname, passwd)
        self.carrier.init()
    def run(self)
        while not self.carrier.done():
            image_name = self.carrier.get_image()
            self.carrier.post_result(self.detect(image_name))
            print "Another image finished!"
        print "All images finished!"

if __name__ == "__main__":
    detector = http_rcnn_detector(prototxt_name, model_name, usrname = usrname, passwd = passwd)
    detector.run()
    
