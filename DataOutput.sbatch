#!/bin/bash
#SBATCH --job-name=DataOutput            # Название задачи
#SBATCH --output=DataOutput.txt       # Файл для вывода результатов
#SBATCH --time=1:00:00                      # Максимальное время выполнения
#SBATCH --cpus-per-task=4                   # Количество CPU на одну задачу
#SBATCH --gpus=0                            # Требуемое кол-во GPU

module load Python/Anaconda_v10.2019        # Загрузка модуля Anaconda
srun python DataOutput.py                # Выполнение расчёта