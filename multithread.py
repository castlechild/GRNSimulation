import multiprocessing
import os


def worker(tasks):
    for task in tasks:
        # Traiter chaque tâche (ici simplement imprimer la tâche)
        print(f"Processing task: {task} sur le coeur{os.getpid()}")


def chunk_tasks(tasks, chunk_size):
    """Divise les tâches en morceaux de taille `chunk_size`."""
    for i in range(0, len(tasks), chunk_size):
        yield tasks[i:i + chunk_size]


if __name__ == "__main__":
    # Création d'une liste de 1000 tâches
    tasks = [f"Task_{i}" for i in range(1000)]

    # Nombre maximum de processus
    num_processes = os.cpu_count()

    # Nombre maximum de tâches par processus
    max_tasks_per_process = 10

    # Diviser les tâches en morceaux
    chunks = list(chunk_tasks(tasks, max_tasks_per_process))

    # Créer et démarrer les processus
    processes = []
    # On ne crée que jusqu'à num_processes processus
    for chunk in chunks:
        # Limiter le nombre de processus en fonction
        # du nombre de cœurs disponibles
        if len(processes) >= num_processes:
            # Attendre que le premier processus se termine
            # avant de créer un nouveau
            processes[0].join()
            processes.pop(0)

        p = multiprocessing.Process(target=worker, args=(chunk,))
        processes.append(p)
        p.start()
    print(os.cpu_count())
    # Attendre la fin de tous les processus
    for p in processes:
        p.join()

    print("All tasks have been processed.")
