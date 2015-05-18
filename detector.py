import edgebox
import nms
import sys
sys.path.insert('/home/maohz12/caffe_test/caffe/python')
import caffe
import httpapi
from PIL import Image
import numpy as np

def rcnn_detect():
    batch_size = 10
    max_num_proposals = 200
    OUTPUT_DIM = 201
    usrname = ''
    passwd = ''
    iou_thres = 0.3

    mean = [104] * 224*224 + [117]*224*224 + [123]*224*224
    mean = np.reshape(mean, (3,224,224), dtype = 'float32')
    net=caffe.Classifier(prototxt_name,model_name,
            mean=mean,
            channel_swap=(2,1,0),
            raw_scale=223,
            image_dims=(224,224))
    caffe.set_mode_gpu()
    httpapi.init(usrname=usrname,passwd=passwd)

    while not httpapi.done():
        image_name = httpapi.get_image()
        image = np.array(Image.load(image_name),dtype='float32')/255
        bboxs = edgebox.get_proposals(image_name, max_num_proposals)
        regions = get_regions(image,bboxs)
        probs = np.zeros(regions.shape[0], OUTPUT_DIM, dtype = 'float32')

        for start in range(0, regions.shape[0], batch_size):
            inputs = np.asarray(map(lambda x:net.transformer.preprocess('data',x),
                regions[start:start+batch_size])) 
            net.forward_all(data=inputs)
            np.copyto(probs[start:start+batch_size], net.blobs['prob'].data)

        keep = nms.nms(bboxs, probs[:,:-1], iou_thres) # Other parameters?
        httpapi.post_result(keep)

if __name__ == "__main__":
    rcnn_detect()
