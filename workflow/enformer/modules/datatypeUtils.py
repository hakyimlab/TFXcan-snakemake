
# Author: Temi
# DAte: Wed Sept 25 2024
# Description: This script is used for ENFORMER inference
# class definitions to pass parameters


import dataclasses
from dataclasses import dataclass, fields
from typing import Union, Any
from datetime import date

# https://stackoverflow.com/questions/56665298/how-to-apply-default-value-to-python-dataclass-field-when-none-was-passed
@dataclass(frozen=True)
class DefaultVal:
    val: Any

@dataclass
class DirectivesHolder:
    # must supply
    model_path: str
    fasta_file: str
    mainRunScript: str
    use_parsl: float
    batch_regions: list
    samples: list
    path_to_vcf: str
    batch_save: bool = True
    prediction_logfile: str = None
    invalid_queries: str = None

    aggregate: Union[None, bool] = None
    aggregate_by_width: bool = True
    aggregation_width: int = None
    aggregation_function: str = None
    bins_indices: Union[int, list] = None
    tracks_indices: Union[int, list] = None
    predictions_expected_shape: Union[tuple, list] = None
    
    # have defaults
    output_file_prefix: str = None
    batch_number: int = 0
    script_path: str = None
    output_directory: str = None
    prediction_logfiles_folder: str = None
    sequence_source: str = None
    debugging: bool = False
    grow_memory: bool = True
    write_log: dict = None
    reverse_complement: bool = False
    tmp_config_path: str = None
    aggregate_by_width: bool = True
    
    write_logdir: str = 'logs'
    run_date: str = DefaultVal(date.today().strftime("%Y-%m-%d"))
    
    

    def __post_init__(self):
        # Loop through the fields
        for field in fields(self):
            # If there is a default and the value of the field is none we can assign a value
            if not isinstance(field.default, dataclasses._MISSING_TYPE) and getattr(self, field.name) is None:
                setattr(self, field.name, field.default)

@dataclass
class ModuleHolder:
    path_to_modules: str
    predictionUtils: str
    checksUtils: str