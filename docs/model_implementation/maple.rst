.. _maple:

Maple
==============

At :code:`tools/ModelGeneration/Maple` we provide the maple Worksheet :code:`CreateModel.mw` which you can use to implement your model and get the tensors.

To add a new model, you have to modify/create five files (for further details, also consult the manual):

1. Go to :code:`include/BSMPT/models` and copy :code:`ClassTemplate.h` to :code:`YourModel.h`. Adjust the name of the class :code:`Class_Template` to :code:`Class_YourModel`.

2. Go to :code:`src/models` and copy :code:`ClassTemplate.cpp` to :code:`YourModel.cpp`, and again change :code:`Class_Template` to :code:`Class_YourModel`. Also, follow the instructions in this file and in the manual to set up your new model. 

3. For your model to compile, you have to open :code:`src/models/CMakeLists.txt` and add :code:`${header_path}/YourModel.h` as well as :code:`YourModel.cpp` to the listed headers and source files.

4. In :code:`include/BSMPT/models/IncludeAllModels.h` you need to add a new entry in the :code:`enum class ModelIDs` above the :code:`stop` entry which is different from the already defined :code:`ModelIDs`, e.g. :code:`YourModel`. Additionally, you have to create a new entry in the :code:`const std::unordered_map<std::string, ModelIDs> ModelNames` map in the same file and add a new line with :code:`{"YourModelName",ModelIDs::YourModel}`.

5. In :code:`src/models/IncludeAllModels.cpp` you have to add :code:`#include <BSMPT/models/YourModel.h>` to the include list. Also, to be able to call your model, you have to extend the :code:`FChoose` function. For this you add a new case to the switch statement, which reads

        case ModelIDs::YourModel: return std::make_unique<Class_YourModel>(); break;

Also contact us if you have a custom model for BSMPT v1.x and you have trouble converting it to the new notation.