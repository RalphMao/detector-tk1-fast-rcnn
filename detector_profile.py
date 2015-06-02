'''
Author: Ralph Mao
May, 2015
'''
import region_proposer
import nms
import cnn
import httpapi
from PIL import Image
import numpy as np
import scipy.misc
import time

# bounding boxes are a list of tuples (xmin, ymin, xmax, ymax)
def get_regions(image, bboxs):                                                                                                        
    image_n = np.array(image)                                                                                                         
    regions = np.zeros((len(bboxs), 224,224,3))                                                                                       
    for i in range(len(bboxs)):                                                                                                       
        regions[i] = scipy.misc.imresize(image_n[bboxs[i,1]:bboxs[i,3], bboxs[i,0]:bboxs[i,2], :], (224 ,224), interp='bilinear')
    return regions                                                                                                                    

class rcnn_detector(object):
    def __init__(self, prototxt_name, model_name, batch_size = 1, max_num_proposals = 60, output_dim = 201, iou_thres = 0.3):
        self.Net = cnn.Net(prototxt_name, model_name, batch_size)
        self.proposer = region_proposer.proposer(max_num_proposals)
        self.nms = nms.reducer(iou_thres)
        self.output_dim = output_dim


    def detect(self, image_in):
        global proposal_time
        global abc_time
        global cnn_time
        global nms_time
        time1 = time.time()
        bboxs = self.proposer.get_proposals(image_in)
        time2 = time.time()
        proposal_time += time2 - time1
        image = Image.open(image_in)
        regions = get_regions(image,bboxs) / 255
        if regions.ndim == 3: # In case of grayscale images
            regions = np.repeat(np.expand_dims(regions,axis=3),3,axis=3)
        probs = np.zeros((regions.shape[0], self.output_dim), dtype = 'float32')
        time3 = time.time()
        abc_time += time3 - time2

        batch_size = self.Net.max_batch
        for start in range(0, regions.shape[0], batch_size):
            np.copyto(probs[start:start + batch_size], self.Net.forward(regions[start:start + batch_size]))
        time4 = time.time()
        cnn_time += time4 - time3

        res = self.nms.multi_class_reduce(bboxs, probs[:,1:]) 
        time5 = time.time()
        nms_time += time5 - time4
        return res



class http_rcnn_detector(rcnn_detector):
    def __init__(self, prototxt_name, model_name, usrname, passwd, batch_size = 1, max_num_proposals = 200, output_dim = 201, iou_thres = 0.3):
        rcnn_detector.__init__(self,prototxt_name, model_name, batch_size, max_num_proposals, output_dim, iou_thres)
        self.carrier = httpapi.carrier(usrname, passwd)
    def run(self):
        while not self.carrier.done():
            image_name = self.carrier.get_image()
            results = self.detect(image_name)
            self.carrier.post_result(results[0],results[1],results[2])
            print "One another image finished!"
        print "All images finished!"

def test_http():
    #==========test http detector============
    usrname = 'nicsefc'
    passwd = 'nics.info'
    detector = http_rcnn_detector(prototxt_name, model_name, batch_size = 10, usrname = usrname, passwd = passwd)
    import time
    start_time = time.time()
    detector.run()
    print "time cost:", time.time() - start_time

def test_mAP():
    #===========test mAP===================
    import time,sys
    detector = rcnn_detector(prototxt_name, model_name, batch_size = 64)
    lines = open('/home/maohz12/DATA/ILSVRC2013_devkit/data/det_lists/val.txt').readlines()
    fout = open('bboxs_results.txt','wb')
    num = 0
    start_time = time.time()
    for line in lines:
        sys.stdout.flush()
        image_name = '/home/maohz12/DATA/ILSVRC2013_DET_val/' + line.split()[0] + '.JPEG'
        image_id = int(line.split()[1])
        class_ids, confidences, bboxs = detector.detect(image_name)
        for i in range(len(class_ids)):
            fout.write('%d %d %f %f %f %f %f\n'%(image_id, class_ids[i], confidences[i], bboxs[i][0], bboxs[i][1], bboxs[i][2], bboxs[i][3]))
        num += 1
        print "Image id%d"%num
    print "Average time cost:%f seconds"%((time.time()-start_time)/num)

if __name__ == "__main__":
    model_name = '/home/maohz12/caffe-tk1/model_test/model/smallVGG11_201/VGG_11_Layer_iter_80000.caffemodel'
    prototxt_name = '/home/maohz12/caffe-tk1/model_test/model/smallVGG11_201/VGG_11_Layer_deploy.prototxt'
    global proposal_time
    global abc_time
    global cnn_time
    global nms_time
    proposal_time = 0.0
    abc_time = 0.0
    cnn_time = 0.0
    nms_time = 0.0
    test_http()
    print proposal_time
    print abc_time
    print cnn_time
    print nms_time

