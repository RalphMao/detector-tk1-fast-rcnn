
import sys
sys.path.insert(0, 'fast-rcnn/lib')
sys.path.insert(0, 'fast-rcnn/caffe-fast-rcnn/python')
import caffe
import fast_rcnn.test
import numpy as np
class fast_rcnn_net(object):
    def __init__(self, prototxt_name, model_name, batch_size = 10, device_id = 0):
        caffe.set_mode_gpu()
        caffe.set_device(device_id)
        self.net = caffe.Net(prototxt_name, model_name, caffe.TEST)

    def detect(self, im, obj_proposals):
        if im == [] or len(obj_proposals) == 0:
            return ([], [])
        return fast_rcnn.test.im_detect(self.net, im, obj_proposals)
        
def test():
    model_name = 'fast-rcnn-model/ilsvrc_fast_rcnn_ft_iter_40000.caffemodel'
    prototxt_name = 'fast-rcnn-model/fast_rcnn_test_new.prototxt'
    net = fast_rcnn_net(prototxt_name, model_name)
    print "Success"

if __name__ == "__main__":
    test()
