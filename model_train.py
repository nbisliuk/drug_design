import keras
from keras import layers
from keras.initializers import RandomNormal

class LSTMModel:

    def __init__(self, model=None):
        self.model = model

    def build_model(self):
        weight_init = RandomNormal(mean=0.0,
                                   stddev=0.05,
                                   seed=11)

        model = keras.models.Sequential()
        model.add(layers.LSTM(256, input_shape=(None, len(word_index) + 1),
                              return_sequences=True, dropout=0.3))

        model.add(layers.LSTM(256, input_shape=(None, len(word_index) + 1),
                              return_sequences=True, dropout=0.5))

        model.add(layers.Dense(len(word_index) + 1, activation='softmax'))

        # optimizer = keras.optimizers.Adam(lr=0.01)
        model.compile(loss='categorical_crossentropy', optimizer='adam')






