
import client
import numpy as np

class carrier(object):
    def __init__(self, usrname, passwd):
        [self.token, self.status] = client.get_token(usrname, passwd)
        [self.total_num, status] = client.get_no_of_images(self.token)
        self.get_idx = 0
        self.post_idx = 0
        print "Total num:", self.total_num
        self.total_num = int(self.total_num)

    def get_image(self):
        if self.done():
            return None
        status = 0
        try:
            status = client.get_image(self.token, self.get_idx + 1, 'images')
            res = 'images/%d.jpg'%(self.get_idx + 1)
        except Exception as e:
            print "Exception:", e
            res = ''
        self.get_idx += 1
        print "Image download. Current idx:", self.get_idx
        return (self.get_idx, res)

    def post_result(self,image_id, class_ids, confidences, bboxs):
        if len(class_ids) != len(confidences) and len(class_ids) != len(bboxs):
            print "Index does not match!"
            return 0

        num_bboxs = len(class_ids)
        if num_bboxs == 0:
            return 1

        data_to_post = {'image_name':[str(image_id)] * num_bboxs,
        'confidence': map(lambda x:str(x), confidences),
        'CLASS_ID': map(lambda x: str(x), class_ids),
        'xmin': map(lambda x: str(x[0]), bboxs),
        'ymin': map(lambda x: str(x[1]), bboxs),
        'xmax': map(lambda x: str(x[2]), bboxs),
        'ymax': map(lambda x: str(x[3]), bboxs)}

        status = 0
        try:
            status = client.post_result(self.token, data_to_post)
        except Exception as e:
            print "Exception:", e
        self.post_idx += 1
        print "Result uploaded. Current idx:", image_id
        return 1

    def done(self):
        if self.status != 1:
            return 1
        return self.total_num == self.get_idx

    def catch(self):
        if self.status != 1:
            return 1
        return self.post_idx == self.get_idx

if __name__ == "__main__":
    import PIL.Image as Image
    import matplotlib.pyplot as plt
    tester = carrier('nicsefc','nics.info')
    print tester.token
    print tester.total_num
    for i in range(15):
        if tester.done():
            break
        image_f = tester.get_image()
        print image_f, i
    ima = Image.open(image_f)
    plt.show(ima)
    print "Get images done"

    bboxs = [(0,0,4,5)]
    confidences = np.array([0.3])
    class_ids = np.array([4])
    for i in range(15):
        tester.post_result(class_ids, confidences, bboxs)

    print "Post result done"
