import numpy as np

class reducer(object):
    def __init__(self, iou_thresh = 0.3, confidence_thresh = 0.00,max_boxs = 5):
        self.iou_thresh = iou_thresh
        self.max_boxs = max_boxs
        self.confidence_thresh = confidence_thresh

    def multi_class_reduce(self, bboxs, probs):
        if len(bboxs) != len(probs):
            return None
        total_num = len(bboxs)
        if total_num == 0:
            return []

        if bboxs.dtype != 'float32':
            bboxs = np.array(bboxs, dtype = 'float32')

        class_id = []
        confidence = np.array([])
        bboxs_chosen = np.array([])
        for i in range(probs.shape[1]):
            scores = probs[:,i]
            keep = self.single_class_reduce(bboxs,scores)
            class_id += [i+1] * len(keep)
            confidence = np.append(confidence, scores[keep])
            bboxs_chosen = np.append(bboxs_chosen, bboxs[keep,:])

        bboxs_chosen = np.reshape(bboxs_chosen, (len(bboxs_chosen)/4, 4))

        class_id = np.array(class_id)
        # Reduce the total number
        I = np.argsort(-confidence)
        return (class_id[I[:5]], confidence[I[:5]], bboxs_chosen[I[:5],:])

    def single_class_reduce(self, bboxs, scores):
        if len(bboxs) != len(scores):
            return None
        total_num = len(bboxs)
        if total_num == 0:
            return []

        pick = []
        idx = np.argsort(-scores)
        idx = idx[(scores[idx] >= self.confidence_thresh)]

        x1 = bboxs[:,0]
        y1 = bboxs[:,1]
        x2 = bboxs[:,2]
        y2 = bboxs[:,3]
        area = (x2 - x1 + 1) * (y2 - y1 + 1)
        while len(idx) > 0 and len(pick) <= self.max_boxs:
            idx_t = idx[0]
            pick.append(idx_t)

            xx1 = np.maximum(x1[idx_t], x1[idx])
            yy1 = np.maximum(y1[idx_t], y1[idx])
            xx2 = np.minimum(x2[idx_t], x2[idx])
            yy2 = np.minimum(y2[idx_t], y2[idx])

            overlap_widths = np.maximum(0, xx2 - xx1 )
            overlap_heights = np.maximum(0, yy2 - yy1 )
            inters = overlap_widths * overlap_heights

            overlap_ratios = inters / (area[idx_t] + area[idx] - inters)
            idx = idx[(overlap_ratios < self.iou_thresh)]

        return pick

if __name__ == "__main__":
    import random
    bboxs = [0,0,20,30] * 5 + [0,0,200,300] * 5 + [40, 50, 250,100] * 5
    bboxs = np.reshape(np.array(bboxs,dtype = 'float32'), (15,4))
    # probs = np.random.uniform(0,1,15*50)
    probs = np.ones((15,50))
    probs = probs.reshape(15,50)
    ducer = reducer()
    a,b,c = ducer.multi_class_reduce(bboxs, probs)
    print a
    print b
    print c
