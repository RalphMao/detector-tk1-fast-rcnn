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
            assert im.shape[1] < 1000
            assert im.shape[2] < 1000
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
        front_process = Process(target=self._front_end, args=(self,0))
        front_process.start()
        front_process2 = Process(target=self._front_end, args=(self,100))
        front_process2.start()

        print "Front process started"
        image_id, im, obj_proposals = self.queue.get()
        start_time = time.time()
        num = 0
        while im != 'Finished':
            self._back_end(image_id, im, obj_proposals)
            image_id, im, obj_proposals = self.queue.get()
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

    def _back_end(self,image_id, im, obj_proposals):
        if im == [] or obj_proposals == []:
            results = ([], [], [])
        else:
            scores, boxes = self.rcnn_net.detect(im, obj_proposals)
            print "Produce %d boxes"%len(scores)
            results = self.reducer.multi_class_reduce(boxes[:,4:], scores[:,1:])
        self.carrier.post_result(image_id, results[0],results[1],results[2])

    @staticmethod
    def _front_end(self,start):
        self.carrier.get_idx = start
        while not self.carrier.done():
            image_id, image_name = self.carrier.get_image()
            print image_name
            try:
                assert os.stat(image_name).st_size < 300 * 1024
                im = cv2.imread(image_name)
                assert len(im.shape) == 3
                assert im.shape[1] < 800
                assert im.shape[2] < 800
            except Exception as e:
                print "Exception:", e
                continue

            obj_proposals = self.proposer.get_proposals(image_name, dtype = 'uint16')
            print "Find %d proposals"%len(obj_proposals)
            if len(obj_proposals) > 0:
                self.queue.put((image_id, im, obj_proposals))

        self.queue.put((0, 'Finish',[]))

def http():
    #==========test http detector for ttq============
    usrname = 'nicsefc'
    passwd = 'nics.info'
    detector = pipe_detector(prototxt_name, model_name, batch_size = 10, usrname = usrname, passwd = passwd)
    detector.run()

if __name__ == "__main__":
    model_name = 'fast-rcnn-model/fast-rcnn-EB-274'
    prototxt_name = 'fast-rcnn-model/test.prototxt'
    '''
    model_name = 'fast-rcnn-model/ilsvrc_fast_rcnn_ft_iter_40000.caffemodel'
    prototxt_name = 'fast-rcnn-model/fast_rcnn_test_new.prototxt'
    '''

    import time, sys
    if len(sys.argv) > 1:
        rcnn = fast_rcnn_detector(prototxt_name, model_name)
        image_in = sys.argv[1]
        results = rcnn.detect(image_in)
        fout = open('res.txt','w')
        for i in range(len(results[0])):
            fout.write('%d %d %f %f %f %f %f\n'%(1, results[0][i], results[1][i], results[2][i][0], results[2][i][1], results[2][i][2], results[2][i][3]))
        fout.close()
    else:
        http()
