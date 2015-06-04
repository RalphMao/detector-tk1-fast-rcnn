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
import os
from multiprocessing import Process, Queue
import time

class fast_rcnn_detector(object):
    def __init__(self, prototxt_name, model_name, batch_size = 1, max_num_proposals = 200, output_dim = 201, iou_thres = 0.3):
        self.proposer = region_proposer.proposer(max_num_proposals)
        self.rcnn_net = fast_rcnn_net.fast_rcnn_net(prototxt_name, model_name, batch_size)
        self.reducer = cnms.reducer(iou_thres)

    def detect(self, image_in):
        try:
            assert os.stat(image_in).st_size < 300 * 1024
            im = cv2.imread(image_in)
            assert len(im.shape) == 3
            assert im.shape[1] < 800
            assert im.shape[2] < 800
        except Exception as e:
            print "Exception:", e
            return ([],[],[])

        obj_proposals = self.proposer.get_proposals(image_in, dtype = 'uint16')
        print "Find %d proposals"%len(obj_proposals)
        scores, boxes = self.rcnn_net.detect(im, obj_proposals)
        print "Produce %d boxes"%len(scores)
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
        print "All images finished!"

class pipe_detector(object):
    def __init__(self, prototxt_name, model_name, usrname, passwd, batch_size = 1, max_num_proposals = 200, output_dim = 201, iou_thres = 0.3):
        self.proposer = region_proposer.proposer(max_num_proposals)
        self.rcnn_net = fast_rcnn_net.fast_rcnn_net(prototxt_name, model_name, batch_size)
        self.reducer = cnms.reducer(iou_thres)
        self.carrier = httpapi.carrier(usrname, passwd)
        self.queue = Queue(8)

        # self._front_end() # Perform first preprocessing for pipeline

    def run(self):
        front_process = Process(target=self._front_end, args=(self,))
        front_process.start()
        print "Front process started"
        im, obj_proposals = self.queue.get()
        start_time = time.time()
        num = 0
        while im != 'Finished':
            self._back_end(im, obj_proposals)
            im, obj_proposals = self.queue.get()
            num += 1
            print "Average speed: %f seconds per image"%((time.time() - start_time) / num)
        print "All images finished"

        front_process.join()


        # im, obj_proposals = self.queue.get()
        # start_time = time.time()
        # num = 0
        # while im != 'Finished':
        #     front_process = Process(self._front_end())
        #     front_process.start()
        #     self._back_end(im, obj_proposals)
        #     im, obj_proposals = self.queue.get()
        #     num += 1
        #     print "Average speed: %f seconds per image"%((time.time() - start_time) / num)
        # print "All images finished"

    def _back_end(self, im, obj_proposals):
        if im == [] or obj_proposals == []:
            results = ([], [], [])
        else:
            scores, boxes = self.rcnn_net.detect(im, obj_proposals)
            print "Produce %d boxes"%len(scores)
            results = self.reducer.multi_class_reduce(boxes[:,4:], scores[:,1:])
            self.carrier.post_result(results[0],results[1],results[2])

    @staticmethod
    def _front_end(self):
        while not self.carrier.done():
            image_name = self.carrier.get_image()
            try:
                assert os.stat(image_name).st_size < 300 * 1024
                im = cv2.imread(image_name)
                assert len(im.shape) == 3
                assert im.shape[1] < 800
                assert im.shape[2] < 800
            except Exception as e:
                print "Exception:", e
                self.queue.put(([],[]))

            obj_proposals = self.proposer.get_proposals(image_name, dtype = 'uint16')
            print "Find %d proposals"%len(obj_proposals)
            self.queue.put((im, obj_proposals))

        self.queue.put(('Finish',[]))




def http():
    #==========test http detector============
    usrname = 'nicsefc'
    passwd = 'nics.info'
    detector = pipe_detector(prototxt_name, model_name, batch_size = 10, usrname = usrname, passwd = passwd)
    detector.run()
def test():
    #===========test mAP===================
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
    model_name = 'fast-rcnn-model/ilsvrc_fast_rcnn_EB_pp_pp_v2_iter_20000.caffemodel'
    prototxt_name = 'fast-rcnn-model/test.prototxt'
    import time, sys
    if len(sys.argv) > 1:
        rcnn = fast_rcnn_detector(prototxt_name, model_name)
        image_in = sys.argv[1]
        results = rcnn.detect(image_in)
        print len(results[0])
    else:
        http()
