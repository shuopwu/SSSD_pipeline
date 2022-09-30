import pickle


def open_custom_input(path):
    with open(path, 'rb') as handle:
            custom_input = pickle.loads(handle.read())
    return custom_input


def save_custom_input(custom_input, save_path):
    with open(save_path,'wb') as outfile:
            pickle.dump(custom_input, outfile,protocol=2)