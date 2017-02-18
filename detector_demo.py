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
from scipy.ndimage import imread

import PIL.Image as Image                                            
import PIL.ImageDraw as Draw                                         
from PIL import ImageFont                                            
def draw(image, bbox, label, size=20):                                        
    draw = Draw.Draw(image)                                          
    draw.line([(bbox[0],bbox[1]), (bbox[0],bbox[3])], fill=(255,0,0))
    draw.line([(bbox[2],bbox[1]), (bbox[2],bbox[3])], fill=(255,0,0))
    draw.line([(bbox[0],bbox[1]), (bbox[2],bbox[1])], fill=(255,0,0))
    draw.line([(bbox[0],bbox[3]), (bbox[2],bbox[3])], fill=(255,0,0))
    draw.rectangle([(bbox[0], bbox[1]), (bbox[0] + int(0.6 * size * len(label)), bbox[1] + size)], fill=(255,255,255))
    font = ImageFont.truetype('font/consolab.ttf', size)               
    draw.text((bbox[0],bbox[1]), label, font=font,fill=(255,0,0))    

def get_synwords(idx):                                               
    file_t = open('font/det_synset_words.txt')                       
    lines = file_t.readlines()                                       
    line = lines[idx]                                                
    word = line.split()[1]                                           
    return word                                                      

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
        front_process = Process(target=self._front_end, args=(self,))
        front_process.start()
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
    def _front_end(self):
        while not self.carrier.done():
            image_id, image_name = self.carrier.get_image()
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
        times = int(sys.argv[1])
    else:
        times = 1
    rcnn = fast_rcnn_detector(prototxt_name, model_name)
    for time_t in range(times, times+40):
        image_in = 'images/%d.jpg'%(time_t+1)
        # add code
        #img_in = cv2.imread("/images/41.jpg")
        img_in = imread('./images/'+str(time_t+1)+'.jpg')
        if len(img_in.shape) != 3: continue
        img_in = img_in[:,:,(2, 1, 0)];
        #img_in = cv2.imread('./images/'+str(time_t+1)+'.jpg')
        #cv2.namedWindow("Image_in")
        #cv2.moveWindow("Image_in",10,0)
        #cv2.imshow("Image_in", img_in)
        #cv2.waitKey(500)
        #cv2.destroyAllWindows()
        # end 
        image = Image.open(image_in)
        results = rcnn.detect(image_in)
        num = 0
        for i in range(len(results[0])):
            word = get_synwords(results[0][i]-1)
            print word
            bbox = results[2][i]
            #if bbox ==  0
            #    break
            draw(image, bbox, word)
            num +=1
            if results[1][i] < 0.25 or num > 3:
                break
       # image.save('bbox_image%d.jpg'%time_t)
    
        image.save('./result/bbox_image%d.jpg'%(time_t+1))
        #img = cv2.imread("./result/bbox_image40.jpg")
        img_out = cv2.imread('./result/bbox_image'+str(time_t+1)+'.jpg')
        cv2.namedWindow("Image_out")
        cv2.moveWindow("Image_out",1000,0)
        cv2.imshow("Image_out", img_out)
        cv2.waitKey(0)
        cv2.destroyAllWindows()
        # image.show()
        print "find %d objects with confidence > 0.25"%num
        # time.sleep(10)
