
import sys
sys.path.insert(0,'/home/maohz12/caffe-tk1/python')
import caffe
import numpy as np

class CNNException(Exception):
    def __init__(self):
        Exception.__init__(self)

def check(condition):
    if not condition:
        raise CNNException

class Net(object):
    def __init__(self, prototxt_name, model_name, max_batch):
        mean = [104] * 224*224 + [117]*224*224 + [123]*224*224
        mean = np.reshape(mean, (3,224,224))
        '''
        lines = open(prototxt_name).readlines()
        lines[2] = "input_dim: %d\n"%max_batch
        f = open('/tmp/.deploy.prototxt','w')
        f.write(''.join(lines))
        '''
        self.net=caffe.Classifier(prototxt_name,model_name,
                mean=mean,
                channel_swap=(2,1,0),
                raw_scale=223,
                image_dims=(224,224))
        self.max_batch = max_batch
        caffe.set_mode_gpu()

    def forward(self, inputs, copy = False):
        check(inputs.shape[0] <= self.max_batch)
        self.net.forward_all(data = np.asarray(map(lambda x:self.net.transformer.preprocess('data',x), inputs)))
        if copy:
            return self.net.blobs['prob'].data[:inputs.shape[0]].copy()
        else:
            return self.net.blobs['prob'].data[:inputs.shape[0]]

def test():
    import PIL.Image as Image
    model_name = '/home/maohz12/ducha-caffe/caffe/smallVGGModel_2_svd_201/snapshot/VGG_11_Layer_iter_80000.caffemodel'
    prototxt_name = '/home/maohz12/ducha-caffe/caffe/smallVGGModel_2_svd_201/VGG_11_Layer_deploy.prototxt'
    net = Net(prototxt_name, model_name, 10)
    print net.net.params['conv1_1'][0].data
    image = Image.open('/home/maohz12/caffe/images/256/ILSVRC2012_val_00002441.jpg')
    image = np.array(image.resize((224,224)))
    inputs = np.array([image,image], dtype='float32')
    print inputs.shape
    probs = net.forward(inputs)
    print probs.shape
    print probs[:,1] 

if __name__ == "__main__":
    test()
