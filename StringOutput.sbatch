#!/bin/bash
#SBATCH --job-name=StringOutput            # Название задачи
#SBATCH --output=StringOutput.txt       # Файл для вывода результатов
#SBATCH --time=6:00:00                      # Максимальное время выполнения
#SBATCH --cpus-per-task=4                   # Количество CPU на одну задачу
#SBATCH --gpus=0                            # Требуемое кол-во GPU

module load Python/Anaconda_v10.2019        # Загрузка модуля Anaconda
srun python StringOutput.py                # Выполнение расчёта