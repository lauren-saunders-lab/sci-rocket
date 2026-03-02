import sys

import functools
import pysam
from frozendict import frozendict
from rapidfuzz import distance, process


# This will allow us to add additional sci-seq related attributes.
class sciRecord:
    def __init__(self, read1: pysam.FastxRecord, read2: pysam.FastxRecord):
        self.read1 = read1
        self.read2 = read2

        # Raw barcode sequences.
        self.p7_sequence, self.p5_sequence = self.read1.comment.split(":")[-1].split("+")
        self.ligation_sequence = None
        self.rt_sequence = None
        self.umi_sequence = None
        self.hashing_sequence = None

        # Names of the barcodes.
        self.p5_name = None
        self.p7_name = None
        self.ligation_name = None
        self.rt_name = None
        self.hashing_name = None

        # Flag of barcode identification.
        self.p5_status = None
        self.p7_status = None
        self.ligation_status = None
        self.rt_status = None
        self.hashing_status = None

        # Sample this read belong to.
        self.sample_name = None
        self.cellular_sequence = None


    def __repr__(self):
        return "@{}:{}+{}\np7:\t{}:{} ({})\np5:\t{}:{} ({})\nlig:\t{}:{} ({})\nrt:\t{}:{} ({})\numi:\t{}\nhash:\t{}:{}({})\nfrom:\t{}".format(
            self.read1.name,
            self.p7_sequence,
            self.p5_sequence,
            self.p7_sequence,
            self.p7_name,
            self.p7_status,
            self.p5_sequence,
            self.p5_name,
            self.p5_status,
            self.ligation_sequence,
            self.ligation_name,
            self.ligation_status,
            self.rt_sequence,
            self.rt_name,
            self.rt_status,
            self.umi_sequence,
            self.hashing_sequence,
            self.hashing_name,
            self.hashing_status,
            self.sample_name,
        )

    # Getter and setters. ------------------------------------------------------
    @property
    def p5_sequence(self):
        return self.__p5_sequence

    @property
    def p7_sequence(self):
        return self.__p7_sequence

    @property
    def ligation_sequence(self):
        return self.__ligation_sequence

    @property
    def rt_sequence(self):
        return self.__rt_sequence

    @property
    def umi_sequence(self):
        return self.__umi_sequence

    @property
    def hashing_sequence(self):
        return self.__hashing_sequence

    @property
    def p5_name(self):
        return self.__p5_name

    @property
    def p7_name(self):
        return self.__p7_name

    @property
    def ligation_name(self):
        return self.__ligation_name

    @property
    def rt_name(self):
        return self.__rt_name

    @property
    def hashing_name(self):
        return self.__hashing_name

    @property
    def p5_status(self):
        return self.__p5_status

    @property
    def p7_status(self):
        return self.__p7_status

    @property
    def ligation_status(self):
        return self.__ligation_status

    @property
    def rt_status(self):
        return self.__rt_status

    @property
    def hashing_status(self):
        return self.__hashing_status

    @property
    def sample_name(self):
        return self.__sample_name

    @property
    def cellular_sequence(self):
        return self.__cellular_sequence

    @p5_sequence.setter
    def p5_sequence(self, p5_sequence):
        self.__p5_sequence = p5_sequence

    @p7_sequence.setter
    def p7_sequence(self, p7_sequence):
        self.__p7_sequence = p7_sequence

    @ligation_sequence.setter
    def ligation_sequence(self, ligation_sequence):
        self.__ligation_sequence = ligation_sequence

    @rt_sequence.setter
    def rt_sequence(self, rt_sequence):
        self.__rt_sequence = rt_sequence

    @umi_sequence.setter
    def umi_sequence(self, umi_sequence):
        self.__umi_sequence = umi_sequence

    @hashing_sequence.setter
    def hashing_sequence(self, hashing_sequence):
        self.__hashing_sequence = hashing_sequence

    @p5_name.setter
    def p5_name(self, p5_name):
        self.__p5_name = p5_name

    @p7_name.setter
    def p7_name(self, p7_name):
        self.__p7_name = p7_name

    @ligation_name.setter
    def ligation_name(self, ligation_name):
        self.__ligation_name = ligation_name

    @rt_name.setter
    def rt_name(self, rt_name):
        self.__rt_name = rt_name

    @hashing_name.setter
    def hashing_name(self, hashing_name):
        self.__hashing_name = hashing_name

    @p5_status.setter
    def p5_status(self, p5_status):
        self.__p5_status = p5_status

    @p7_status.setter
    def p7_status(self, p7_status):
        self.__p7_status = p7_status

    @ligation_status.setter
    def ligation_status(self, ligation_status):
        self.__ligation_status = ligation_status

    @rt_status.setter
    def rt_status(self, rt_status):
        self.__rt_status = rt_status

    @hashing_status.setter
    def hashing_status(self, hashing_status):
        self.__hashing_status = hashing_status

    @sample_name.setter
    def sample_name(self, sample_name):
        self.__sample_name = sample_name

    @cellular_sequence.setter
    def cellular_sequence(self, cellular_sequence):
        self.__cellular_sequence = cellular_sequence

    # Methods -----------------------------------------------------------------

    def determine_p5(self, p5_barcodes):
        """
        Determine the p5 barcode of the read (reverse complement) within the read-name.

        Parameters:
            p5_barcodes (dict): Dictionary of p5 barcodes.

        Sets:
            p5_name (str): Name of the p5 barcode.
            p5_status (str): Status of the p5 barcode.
            p5_sequence (str): Sequence of the p5 barcode (if rescued)
        """
        try:
            self.p5_name = p5_barcodes[self.p5_sequence]
            self.p5_status = "Correct"
        except KeyError:
            p5_sequence, self.p5_name = self.__find_closest_match(self.p5_sequence, p5_barcodes)
            if self.p5_name != None:
                self.p5_status = "Corrected"
                self.p5_sequence = p5_sequence 

    def determine_p7(self, p7_barcodes):
        """
        Determine the p7 barcode of the read within the read-name.

        Parameters:
            p7_barcodes (dict): Dictionary of p7 barcodes.

        Sets:
            p7_name (str): Name of the p7 barcode.
            p7_status (str): Status of the p7 barcode.
            p7_sequence (str): Sequence of the p7 barcode (if rescued)
        """
        try:
            self.p7_name = p7_barcodes[self.p7_sequence]
            self.p7_status = "Correct"
        except KeyError:
            p7_sequence, self.p7_name = self.__find_closest_match(self.p7_sequence, p7_barcodes)
            if self.p7_name != None:
                self.p7_status = "Corrected"
                self.p7_sequence = p7_sequence

    def determine_ligation(self, ligation_barcodes):
        """
        Determine the ligation barcode of the read within R1.
        This ligation barcode varies in length from 9 to 10nt.. So check for both.

        Parameters:
            ligation_barcodes (dict): Dictionary of ligation barcodes.

        Sets:
            ligation_name (str): Name of the ligation barcode.
            ligation_status (str): Status of the ligation barcode.
            ligation_sequence (str): Sequence of the ligation barcode.
        """
        sequence_ligation_9nt = self.read1.sequence[0:9]
        sequence_ligation_10nt = self.read1.sequence[0:10]
        self.ligation_sequence = sequence_ligation_10nt

        # First try exact match 10nt.
        try:
            self.ligation_name = ligation_barcodes["ligation_10nt"][sequence_ligation_10nt]
            self.ligation_status = "Correct"
            self.ligation_sequence = sequence_ligation_10nt
        except KeyError:
            # If no exact match, try 9nt.
            try:
                self.ligation_name = ligation_barcodes["ligation_9nt"][sequence_ligation_9nt]
                self.ligation_status = "Correct"
                self.ligation_sequence = sequence_ligation_9nt
            except KeyError:
                # If no exact match, try hamming distance of 1.
                ligation_sequence, self.ligation_name = self.__find_closest_match(sequence_ligation_10nt, ligation_barcodes["ligation_10nt"])
                if self.ligation_name != None:
                    self.ligation_status = "Corrected"
                    self.ligation_sequence = ligation_sequence
                else:
                    ligation_sequence, self.ligation_name = self.__find_closest_match(sequence_ligation_9nt, ligation_barcodes["ligation_9nt"])
                    if self.ligation_name != None:
                        self.ligation_status = "Corrected"
                        self.ligation_sequence = ligation_sequence

    def determine_rt(self, rt_barcodes):
        """
        Determines the RT barcode of the read within R1.
        This RT barcode is always 10nt long, based on the ligation it can be the last 10nt or the last 10nt minus one.

        Parameters:
            rt_barcodes (dict): Dictionary of RT barcodes.

        Sets:
            rt_name (str): Name of the RT barcode.
            rt_status (str): Status of the RT barcode.
            rt_sequence (str): Sequence of the RT barcode.
        """
        # Check length of the ligation barcode to determine the location of the other barcodes.
        if self.ligation_status == None or len(self.ligation_sequence) == 10:
            # Retrieve the RT barcode from R1 (last 10 bp).
            self.rt_sequence = self.read1.sequence[24:34]
        else:
            # Retrieve the RT barcode from R1 (last 10 bp, minus one).
            self.rt_sequence = self.read1.sequence[23:33]
        
        try:
            self.rt_name = rt_barcodes[self.rt_sequence]
            self.rt_status = "Correct"
        except KeyError:
            rt_sequence, self.rt_name = self.__find_closest_match(self.rt_sequence, rt_barcodes)
            if self.rt_name != None:
                self.rt_status = "Corrected"
                self.rt_sequence = rt_sequence

    def determine_umi(self):
        """
        Determine the UMI of the read within R1.
        This UMI is always 8nt long and is located after the primer sequence of 6nt, based on the length of the ligation barcode.

        Sets:
            umi_sequence (str): Sequence of the UMI.
        """
        # Check length of the ligation barcode to determine the location of the other barcodes.
        if self.ligation_status == None or len(self.ligation_sequence) == 10:
            # Retrieve the UMI from R1 (next 8 bp after the primer sequence of 6nt).
            self.umi_sequence = self.read1.sequence[16:24]
        else:
            # Retrieve the UMI from R1 (next 8 bp after the primer sequence of 6nt).
            self.umi_sequence = self.read1.sequence[15:23]

    def determine_sample(self, dict_samples):
        """
        Determine the sample this read belongs to based on p5, p7 and RT.
        Also sets the cellular barcode based on p5, p7, ligation and RT.

        Parameters:
            dict_samples (dict): Dictionary of samples.

        Sets:
            sample_name (str): Name of the sample.
            cellular_sequence (str): Sequence of the cellular barcode.
        """
        if self.p5_name != None and self.p7_name != None and self.rt_name != None and self.ligation_name != None:
            try:
                self.sample_name = dict_samples[self.p5_name + "_" + self.p7_name + "_" + self.rt_name]
            except KeyError:
                self.sample_name = None

        # Set the cellular sequence (sequence)
        if self.sample_name != None:
            if len(self.ligation_sequence) == 10:
                self.cellular_sequence = self.p7_sequence + self.p5_sequence + self.ligation_sequence + self.rt_sequence
            else:
                self.cellular_sequence = self.p7_sequence + self.p5_sequence + self.ligation_sequence + "G" + self.rt_sequence

    def convert_R1andR2(self):
        """
        Convert the original R1 sequence to the new R1 sequence.
        This new R1 sequence contains the cellular barcode sequences + UMI.

        The read-names of R1 and R2 are also updated to contain the barcodes.

        Sets:
            read1.sequence (str): New R1 sequence.
            read2.sequence (str): New R2 sequence.
            read1.name (str): New R1 name.
            read2.name (str): New R2 name.
        """
        self.read1.sequence = self.cellular_sequence + self.umi_sequence
        self.read1.quality = "F" * len(self.read1.sequence)

        self.read1.set_name("{}|P5{}-P7{}|{}|{}_{}".format(self.read1.name, self.p5_name, self.p7_name, self.ligation_name, self.rt_name, self.umi_sequence))
        self.read2.set_name("{}|P5{}-P7{}|{}|{}_{}".format(self.read2.name, self.p5_name, self.p7_name, self.ligation_name, self.rt_name, self.umi_sequence))

    def determine_hash(self, dict_hashing):
        """
        If the sample is known to contain hash-barcodes:
        - Check for the presence of a polyA.
        - Next, check if the 10nt hash + rescue (1 hamming) is directly upstream of polyA (-1nt spacer).
        - If not, check for a perfect match of any hash-barcode in R2 (pre-compiled regex) prior to polyA.
        - After counting, these hash read-pairs are discarded.

        Parameters:
            dict_hashing (dict): Dictionary of hashing barcodes.

        Sets:
            hashing_name (str): Name of the hashing barcode.
            hashing_status (str): Status of the hashing barcode.
            hashing_sequence (str): Sequence of the hashing barcode.
        """
        if self.sample_name != None and dict_hashing and self.sample_name in dict_hashing:
            # Check for polyA in R2.
            start_poly = str.find(self.read2.sequence, "AAAA")

            if start_poly != -1:
                # Check hash-barcode upstream of polyA - 1nt.
                if start_poly - 11 < 0:
                    sequence_hash_raw = self.read2.sequence[0 : start_poly - 1]
                else:
                    sequence_hash_raw = self.read2.sequence[start_poly - 11 : start_poly - 1]

                try:
                    self.hashing_name = dict_hashing[self.sample_name]["sheet"][sequence_hash_raw]
                    self.hashing_sequence = sequence_hash_raw
                    self.hashing_status = "Correct"
                except KeyError:
                    # Rescue the sequence directly prior to polyA.
                    self.hashing_sequence, self.hashing_name = self.__find_closest_match(sequence_hash_raw, dict_hashing[self.sample_name]["sheet"])

                    if self.hashing_name != None:
                        self.hashing_status = "Corrected"
                    else:
                        # Check for the presence of any hash in the entire R2 sequence prior to poly-A.
                        result = dict_hashing[self.sample_name]["regex"].findall(self.read2.sequence[:start_poly])

                        if len(result) == 1:
                            self.hashing_sequence = result[0]
                            self.hashing_name = dict_hashing[self.sample_name]["sheet"][self.hashing_sequence]
                            self.hashing_status = "Corrected upstream"

    # Define hidden methods. ---------------------------------------------------

    @staticmethod
    @functools.lru_cache(maxsize=None)
    def __find_closest_match(sequence: str, comparison: frozendict):
        """
        Find the closest match between a sequence and a dictionary of sequences.
        If there is only one match with a hamming distance of 1, return the name and sequence.
        If there are multiple matches with a hamming distance of 1, return None.

        Parameters:
            sequence (str): Sequence to compare.
            comparison (dict): Dictionary of sequences to compare to.

        Returns:
            str: Sequence of the closest match.
            str: Name of the closest match.
        """

        # Calculate the hamming distance between sequence and all keys in dict.
        # Only keep the keys with a hamming distance of 1.
        distances = process.extract(sequence, comparison.keys(), scorer=distance.Hamming.distance, score_cutoff=1, limit=2)

        # If there is only one key with a distance of 1, return the name and sequence.
        if len(distances) == 1:
            sequence = distances[0][0]
            return sequence, comparison[sequence]
        else:
            return None, None
