
import os
import sys
import pandas as pd


########################################
### file create and search functions ###
########################################


def find_file(path_b,pattern_list = []):
    
    listfiles_path = []
    for path, dirs, files in os.walk(path_b):
        for f in files:
            if all([p in f for p in pattern_list]):
                listfiles_path.append(os.path.join(os.path.realpath(path), f))

    files_with_path = sorted(listfiles_path)

    return files_with_path

def create_folder(path, name):
    filename = os.path.join(path, name)
    if not os.path.isdir(filename):
        os.makedirs(filename)
    
    return filename


def find_files_in_list(f_list,keyword_list=[]):
    retult_f = []
    for f in f_list :
        if all([keyword in f for keyword in keyword_list]):
            retult_f.append(f)
    return sorted(retult_f)


####### count rows #########
def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b

def count_rows(files):
    d = {'file':[], 'how_many_rows':[]}
    for file in files:
        with open(file, "r") as f:
            d['file'].append(os.path.basename(file))

            count = sum(bl.count("\n") for bl in blocks(f)) - 1
            d['how_many_rows'].append(count) 
            
    return pd.DataFrame(d)

######### save dictionary to json #########
import json
class JSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if hasattr(obj, 'to_json'):
            return obj.to_json(orient='records')
        return json.JSONEncoder.default(self, obj)
    
def save_dict_pandas_helper(data,save_dir, save_name):
    with open(os.path.join(save_dir,save_name), 'w') as fp:
        json.dump(data, fp, cls=JSONEncoder)
    print('{}'.format(os.path.join(save_dir,save_name)))
    return os.path.join(save_dir,save_name)

######### xlxs #########

def append_xlxs(data,f,sheet_name):
    from openpyxl import load_workbook
    writer = pd.ExcelWriter(f, engine = 'openpyxl')
    writer.book = load_workbook(f)

    data.to_excel(writer, sheet_name = sheet_name,index=False)
    
    writer.save()
    writer.close()
    print('save {} in {}'.format(sheet_name,f))
    
    
    
if __name__ == " __main__" :
    # below is an example : how to use append xlxs with pandas dataframe
    # read xlxs
    f = 'example.xlxs'
    with pd.ExcelWriter(f,engine='xlsxwriter') as writer:
        data1.to_excel(writer,sheet_name='table_name1',index=False)
    append_xlxs(data2,f,'table_name2')
    
    # read dataframe from f 
    db = pd.ExcelFile(f)
    db.sheet_names
    data = db.parse('table_name')

