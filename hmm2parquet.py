from ast import parse
import pandas as pd

class ProfilesFile:
    def __init__(self):
        self.profiles = {}

    def read_file(self, filepath):
        if not filepath.endswith(".hmm"):
            raise ValueError("File must be a .hmm file")

        with open(filepath, "r") as f:
            curr_profile = None
            header_mode = True
            row_mode = "match"  # one of match, insert, delete
            for line in f:
                line = line.strip()
                if header_mode:
                    if line.startswith("HMMER3"):
                        continue
                    if line.startswith("NAME"):
                        print(f"Processing profile no. {1 + len(self.profiles)}...")
                        name = line.split()[1]
                        curr_profile = ProfileHMM()
                        curr_profile.header["name"] = name
                        continue
                    if line.startswith("HMM"):
                        header_mode = False
                        continue
                    else:
                        key, value = self._parse_key_value(line)
                        curr_profile.header[key.lower()] = value
                else:
                    # not in header mode -> table mode
                    if line == ("//"):
                        self.profiles[curr_profile.header['name']] = curr_profile
                        curr_profile.save_to_file()
                        curr_profile = None
                        header_mode = True
                        continue
                    if line.startswith("m->m"):
                        continue
                    else:
                        if row_mode == "match":
                            curr_profile.insert_row(line, "match")
                            row_mode = "insert"
                        elif row_mode == "insert":
                            curr_profile.insert_row(line, "insert")
                            row_mode = "delete"
                        elif row_mode == "delete":
                            curr_profile.insert_row(line, "delete")
                            row_mode = "match"
                        else:
                            raise ValueError("Invalid row mode") 
                        
    def to_parquet(self, folderpath):
        for name, profile in self.profiles.items():
            profile.to_parquet(folderpath)

                            


    def _parse_key_value(self, str):
        return (str.split()[0].strip(), " ".join(str.split()[1:]).strip())


class ProfileHMM:
    """Represents a single profile HMM"""

    def __init__(
        self,
        header = {}
    ):
        self.header = header
        self.table = pd.DataFrame(columns=["state", "mode", "A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","m->m","m->i","m->d","i->m","i->i","d->m","d->d", "rest"])
        self.curr_state = None
        self.data = []

    def save_to_file(self):
        self.table = self.create_df_from_data()
        self.table.to_parquet(f"pfam_output/{self.header['name']}.parquet")

    def create_df_from_data(self):
        return pd.DataFrame(self.data, columns=["state","mode","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","m->m","m->i","m->d","i->m","i->i","d->m","d->d", "rest"])

    def insert_row(self, line, row_mode:str):
        with_rest = None
        if row_mode == "match":
            with_rest = self._parse_normal_row(line, row_mode)
        elif row_mode == "insert":
            with_rest = self._parse_normal_row(line, row_mode)
        elif row_mode == "delete":
            with_rest = self._parse_delete_row(line, row_mode)
        else:
            raise ValueError("Invalid row mode")     

        #row = pd.DataFrame([with_rest], columns=["state","mode","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","m->m","m->i","m->d","i->m","i->i","d->m","d->d", "rest"])
        #self.table = self.table.append(row)
        self.data.append(with_rest)

    def _parse_normal_row(self, line, row_mode):
        parsed = [el.strip() for el in line.split()]
        if row_mode == "match":
            self.curr_state = parsed[0]
        if row_mode == "insert":
            parsed.insert(0, self.curr_state)
        parsed.insert(1, row_mode)
        for i in range(22):
            if i > 1:
                parsed[i] = float(parsed[i])
        # pop rest into last column
        rest = parsed[22:]            
        rest_str = "".join(rest)
        # add None entries for transition probabilities
        without_rest = parsed[:22]
        for i in range(7):
            without_rest.insert(i+22, None)
        with_rest = without_rest + [rest_str]
        return with_rest
    
    def _parse_delete_row(self, line, row_mode):
        parsed = [el.strip() for el in line.split()]
        total = [self.curr_state, row_mode]

        for i in range(21):
            total.append(None)            
        
        # add None entries for transition probabilities
        for i in range(7):
            try:
                total.insert(i+22, float(parsed[i]))
            except:
                total.insert(i+22, -1)
        
        return total
    
    def to_parquet(self, folderpath):
        self.table.to_parquet(f"{folderpath}/{self.header['name']}.parquet")