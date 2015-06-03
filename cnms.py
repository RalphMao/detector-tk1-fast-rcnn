
import sys
sys.path.insert(0, 'fast-rcnn/lib')
import utils.cython_nms
import numpy as np

class reducer(object):
    def __init__(self, iou_thresh = 0.3, confidence_thresh = 0.01,max_boxs = 40):
        self.iou_thresh = iou_thresh
        self.max_boxs = max_boxs
        self.confidence_thresh = confidence_thresh
        self.max_boxs_per_class = max_boxs / 2 # Heuristic setting

    def multi_class_reduce(self, boxes, scores):
        if boxes.shape[0] != scores.shape[0] or boxes.shape[1] != scores.shape[1] * 4:
            print "boxes shape:", boxes.shape
            print "scores shape:", scores.shape
            raise ValueError('Index does not match!')

        num_boxes = boxes.shape[0]
        num_cls = scores.shape[1]
        class_id = []
        confidence = np.array([], dtype = 'float32')
        bboxs_chosen = np.array([], dtype = 'float32')
        for cls_ind in range(num_cls):
            cls_boxes = boxes[:, 4*cls_ind:4*(cls_ind + 1)]
            cls_scores = scores[:, cls_ind]
            keep_c = np.where(cls_scores > self.confidence_thresh)[0]
            cls_boxes = cls_boxes[keep_c, :]
            cls_scores = cls_scores[keep_c]

            dets = np.hstack((cls_boxes, cls_scores[:, np.newaxis])).astype(np.float32)

            keep = utils.cython_nms.nms(dets, self.iou_thresh)
            '''
            keep2 = np.where(scores[:, cls_ind] > self.confidence_thresh)[0]
            keep = np.intersect1d(keep1, keep2)
            '''
            if len(keep) == 0:
                continue

            class_id += [cls_ind + 1] * len(keep)
            confidence = np.append(confidence, cls_scores[keep])
            bboxs_chosen = np.append(bboxs_chosen, cls_boxes[keep,:])

        bboxs_chosen = np.reshape(bboxs_chosen, (len(bboxs_chosen)/4, 4))
        class_id = np.array(class_id)
        I = np.argsort(-confidence)
        return (class_id[I[:self.max_boxs]], confidence[I[:self.max_boxs]], bboxs_chosen[I[:self.max_boxs],:])

if __name__ == "__main__":
    import random
    bboxs = [0,0,20,30] * 5 + [0,0,200,300] * 5 + [40, 50, 250,100] * 5
    bboxs = np.reshape(np.array(bboxs * 50,dtype = 'float32'), (15,200))
    # probs = np.random.uniform(0,1,15*50)
    probs = np.ones((15,50))
    probs = probs.reshape(15,50)
    ducer = reducer(max_boxs = 5)
    a,b,c = ducer.multi_class_reduce(bboxs, probs)
    print a
    print b
    print c
