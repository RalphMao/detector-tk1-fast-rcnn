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
        self.reducer = cnms.reducer(iou_thres)

    def detect(self, image_in):
        global proposal_time
        global cnn_time
        global nms_time
        time1 = time.time()
        obj_proposals = self.proposer.get_proposals(image_in, dtype = 'uint16')
        im = cv2.imread(image_in)
        time2 = time.time()
        proposal_time += time2 - time1
        scores, boxes = self.rcnn_net.detect(im, obj_proposals)
        time3 = time.time()
        cnn_time += time3 - time2
        ans = self.reducer.multi_class_reduce(boxes[:,4:], scores[:,1:])
        time4 = time.time()
        nms_time += time4 - time3
        return ans


class http_rcnn_detector(fast_rcnn_detector):
    def __init__(self, prototxt_name, model_name, usrname, passwd, batch_size = 1, max_num_proposals = 200, output_dim = 201, iou_thres = 0.3):
        fast_rcnn_detector.__init__(self,prototxt_name, model_name, batch_size, max_num_proposals, output_dim, iou_thres)
        self.carrier = httpapi.carrier(usrname, passwd)
    def run(self):
        global down_time
        global post_time
        max_num = 200
        while not self.carrier.done() and max_num > 0:
            start_time = time.time()
            image_name = self.carrier.get_image()[1]
            down_time += time.time() - start_time
            results = self.detect(image_name)
            start_time = time.time()
            self.carrier.post_result(0,results[0],results[1],results[2])
            post_time += time.time() - start_time
            print "One another image finished!"
            max_num -=1
        print "All images finished!"

def http(num_proposals):
    #==========test http detector============
    usrname = 'nicsefc'
    passwd = 'nics.info'
    detector = http_rcnn_detector(prototxt_name, model_name, batch_size = 10, usrname = usrname, passwd = passwd, max_num_proposals=num_proposals)
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

    
def multi_profile():
    num_list = [20,50,100,150,200,400,700,1000]
    fout = open('frtime_num.log','w')
    global proposal_time
    global cnn_time
    global nms_time
    global down_time
    global post_time
    for num in num_list:
        single_profile(num)
        fout.write('%d, %f\n'%(num, cnn_time))

    fout.close()
def single_profile(num_proposals):
    global proposal_time
    global cnn_time
    global nms_time
    global down_time
    global post_time
    down_time = 0.0
    post_time = 0.0
    proposal_time = 0.0
    cnn_time = 0.0
    nms_time = 0.0
    http(num_proposals)
    print proposal_time
    print cnn_time
    print nms_time
    print down_time
    print post_time

if __name__ == "__main__":
    model_name = 'fast-rcnn-model/fast-rcnn-EB-274'
    prototxt_name = 'fast-rcnn-model/test.prototxt'
    import time

    multi_profile()
