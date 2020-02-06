import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
import pandas as pd
import yaml
import os

with open("config.yml", 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

#import data
bulkdir = "Data/Slideseq/" + cfg['default']['slideseqfolder'] + "/results/Bulk"
X = pd.read_csv(os.path.join(bulkdir, "X_bulk.csv"))
b = pd.read_csv(os.path.join(bulkdir, "b_bulk.csv"))
cell_type_names = X.columns[1:]
X = X.iloc[:,1:].to_numpy()
b = b.iloc[:,1].to_numpy()

#define the model
class Model(object):
  def __init__(self):
    self.N = len(cell_type_names)
    self.w = tf.Variable([1/self.N for i in range(self.N)])

  def __call__(self, x):
    return tf.tensordot(x[0:(self.N)],tf.abs(self.w),1)
    
  def predict_batch(self, x):
    return tf.linalg.matvec(x[:,0:(self.N)],tf.abs(self.w)) 

def loss(predicted_y, target_y):
  return tf.math.reduce_sum(tf.square(tf.math.log(predicted_y) - tf.math.log(target_y)))

model = Model()

#create the loss function
X_var = tf.compat.v1.placeholder("float32") 
Y_var = tf.compat.v1.placeholder("float32") 
my_loss = loss(model(X_var),Y_var)
batch_loss = loss(model.predict_batch(X_var),Y_var)
optimizer = tf.train.AdamOptimizer(1e-4).minimize(my_loss)
init = tf.global_variables_initializer() 


#train the model
training_epochs = 50
with tf.Session() as sess: 
     
    # Initializing the Variables 
    sess.run(init) 
      
    # Iterating through all the epochs 
    for epoch in range(training_epochs): 
          
        # Feeding each data point into the optimizer using Feed Dictionary 
        for (_x, _y) in zip(X, b): 
            sess.run(optimizer, feed_dict = {X_var : _x, Y_var : _y}) 
          
        # Displaying the result after every 5 epochs 
        if (epoch + 1) % 5 == 0: 
            # Calculating the cost a every epoch 
            c = sess.run(batch_loss, feed_dict = {X_var : X, Y_var : b}) 
            print("Epoch", (epoch + 1), ": cost =", c, "w =", sess.run(model.w)) 
      
    # Storing necessary values to be used outside the Session 
    training_cost = sess.run(batch_loss, feed_dict = {X_var : X, Y_var : b}) 
    weight = sess.run(model.w) 

#save the results
weight_df = pd.DataFrame(abs(weight),index=cell_type_names,columns=["Weight"])
weight_df.to_csv(os.path.join(bulkdir, "weights.csv"))
