# coding: utf-8
# exception_with_code.py

from abc import ABC, abstractmethod

import error_codes


#################### ExceptionWithCode ####################

class ExceptionWithCode(Exception, ABC):
	@abstractmethod
	def __init__(self, message: str):
		super().__init__(message)
	
	@abstractmethod
	def get_error_code(self) -> error_codes.Type:
		pass


#################### FileNotFoundException ####################

class FileNotFoundException(ExceptionWithCode):
	def __init__(self, path: str):
		super().__init__("error : " + path + " can't open.")
	
	def get_error_code(self) -> error_codes.Type:
		return error_codes.Type.FILE_NOT_FOUND
