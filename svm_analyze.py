import mne
mne.set_log_level('error')
import numpy as np
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import classification_report
from sklearn.metrics import accuracy_score, precision_score, recall_score
from sklearn.metrics import confusion_matrix
import sys
from itertools import islice\



def separate_sentences(eeg_data, vmrk_file, sentence_num = 639, sentence_length = 1281): # input .vhdr file, desired # of sentences to seperate

    eeg_data = eeg_data._data
    
    seperated_data = [np.empty([eeg_data.shape[0], sentence_length]) for i in range(sentence_num)] # initialize empty training matrices
    i = 0

    with open(vmrk_file) as file:
    
        for j in range(13): # skip the header
            next(file)
            
        for line in file:
            time = line.split(",")[2] # get second number (time marker)
            start = int(time) + 100 # 100 ms offset
            stop = start + sentence_length
            seperated_data[i] = eeg_data[16:, start:stop,]
            i += 1

    return seperated_data # list of np training matrices


def generate_labels(dat_file, sentence_num = 639): # input .dat file, return one hot encoding of data class

    labels = [np.empty(4) for i in range(sentence_num)]
    sentence = 0

    with open(dat_file) as file:

        for i in range(15): # skip the header
            next(file)

        for line in file:
            
            code = int(line.split()[1]) # get second number (the sentence code)

            if code < 13: # balanced
                labels[sentence] = 1

            elif code < 24: # unbalanced
                labels[sentence] = 2

            elif code < 37: # imperative
                labels[sentence] = 3
            
            elif line.split()[2] == "blah" : # blahs
                labels[sentence] = 4
            
            else: # story
                labels[sentence] = 5

            sentence += 1

    return labels


if __name__ == "__main__":

    data = mne.io.read_raw_brainvision("parser994.vhdr", preload = True)

    filt_data = data.filter(0.1, 30, method = "iir") # bandpass filter, iir/butterworth
   
    separated_data = separate_sentences(filt_data, "parser994.vmrk")

    freq_data = []
    for item in separated_data:
        fft_ed = np.fft.fft2(item).flatten(order='F')
        fft_ed_append = np.concatenate((fft_ed.real, fft_ed.imag))
        freq_data.append(fft_ed_append)

    labels = generate_labels("994.dat")

    #train_data, test_data, train_labels, test_labels = train_test_split(freq_data, labels, test_size=0.33, random_state=42)

    svm = SVC(kernel='rbf', C = 10, gamma = 1) # hyper parameters tuned with GridSearch

    # Use cross-validation to generate predictions for each fold
    predicted_labels = cross_val_predict(svm, freq_data, labels, cv=5)

    report = classification_report(labels, predicted_labels)

    print(report)