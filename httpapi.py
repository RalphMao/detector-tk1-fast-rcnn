
import client

class carrier(object):
    def __init__(self, usrname, passwd):
        self.token = client.get_token(usrname, passwd)
        self.total_num = client.get_no_of_images(self.token)
        self.get_idx = 0
        self.post_idx = 0

    def get_image(self):
        if self.done():
            return None
        status = 0
        while status == 0:
            try:
                status = client.get_image(self.token, self.get_idx + 1, 'images')
            except Exception as e:
                print "Exception:", e
                print "Try to fetch the image %d again"%(self.get_idx + 1)
        self.get_idx += 1
        return 'images/%d.jpg'%self.get_idx

    def post_result(self,class_ids, confidences, bboxs):
        if self.catch():
            return 0
        if len(class_ids) != len(confidences) and len(class_ids) != len(bboxs):
            return 0

        num_bboxs = len(class_ids)
        data_to_post = {'image_name':[str(self.post_idx+1)] * num_bboxs,
        'confidence': map(lambda x:str(x), confidences),
        'CLASS_ID': map(lambda x: str(x), class_ids),
        'xmin': map(lambda x: str(x[0]), bboxs),
        'ymin': map(lambda x: str(x[1]), bboxs),
        'xmax': map(lambda x: str(x[2]), bboxs),
        'ymax': map(lambda x: str(x[3]), bboxs)}

        status = 0
        while status == 0:
            try:
                status = client.post_result(self.token, data_to_post)
            except Exception as e:
                print "Exception:", e
                print "Try to fetch the image %d again"%(self.get_idx + 1)
        self.post_idx += 1
        return 1

    def done(self):
        return self.total_num == self.get_idx

    def catch(self):
        return self.post_idx == self.get_idx
