# coding: utf-8
# error_codes.py

from enum import Enum


#################### Type ####################

class Type(Enum):
	NO_ERROR = 0
	WRONG_ARGUMENT = 1
	FILE_NOT_FOUND = 134
	PEDIGREE_INVALID_FORMAT = 135
	PARENT_NOT_DEFINED = 136
	SAMPLES_NOT_IN_PEDIGREE = 137
	VCF_INVALID_FORMAT = 139
	MAP_INVALID_FORMAT = 140
