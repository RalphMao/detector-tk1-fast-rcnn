'''
Author: Ralph Mao
May, 2015
'''
import region_proposer
import httpapi
import fast_rcnn_net
import cnms
import cv2
import numpy as np

class fast_rcnn_detector(object):
    def __init__(self, prototxt_name, model_name, batch_size = 1, max_num_proposals = 200, output_dim = 201, iou_thres = 0.3):
        self.proposer = region_proposer.proposer(max_num_proposals)
        self.rcnn_net = fast_rcnn_net.fast_rcnn_net(prototxt_name, model_name, batch_size)
        self.reducer = cnms.reducer(iou_thres, max_num_proposals)

    def detect(self, image_in):
        obj_proposals = self.proposer.get_proposals(image_in, dtype = 'uint16')
        im = cv2.imread(image_in)
        scores, boxes = self.rcnn_net.detect(im, obj_proposals)
        ans = self.reducer.multi_class_reduce(boxes[:,4:], scores[:,1:])
        return ans


class http_rcnn_detector(fast_rcnn_detector):
    def __init__(self, prototxt_name, model_name, usrname, passwd, batch_size = 1, max_num_proposals = 200, output_dim = 201, iou_thres = 0.3):
        fast_rcnn_detector.__init__(self,prototxt_name, model_name, batch_size, max_num_proposals, output_dim, iou_thres)
        self.carrier = httpapi.carrier(usrname, passwd)
    def run(self):
        while not self.carrier.done():
            image_name = self.carrier.get_image()
            results = self.detect(image_name)
            self.carrier.post_result(results[0],results[1],results[2])
            print "One another image finished!"
        print "All images finished!"

def http():
    #==========test http detector============
    usrname = 'nicsefc'
    passwd = 'nics.info'
    detector = http_rcnn_detector(prototxt_name, model_name, batch_size = 10, usrname = usrname, passwd = passwd)
    detector.run()
def test():
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
    model_name = '/home/maohz12/rcnn-python/fast-rcnn-model/ilsvrc_fast_rcnn_ft_iter_40000.caffemodel'  
    prototxt_name = '/home/maohz12/rcnn-python/fast-rcnn-model/fast_rcnn_test_new.prototxt'
    http()
