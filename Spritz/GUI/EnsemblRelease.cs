using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;

namespace Spritz
{
    public class EnsemblRelease
    {
        public string Release { get; set; }
        public ObservableCollection<string> Species { get; set; }
        public Dictionary<string, string> Genomes { get; set; } // Mus_musculus GRCm38
        public Dictionary<string, string> Organisms { get; set; } // Mus_musculus GRCm38

        public static ObservableCollection<EnsemblRelease> GetReleases()
        {
            var ensemblReleases = new ObservableCollection<EnsemblRelease>();

            // read release.txt files into a list
            var genomeDB = File.ReadAllLines(Path.Combine(Directory.GetCurrentDirectory(), "genomes.csv")).Where(line => !line.StartsWith("#")).ToList();
            var releases = genomeDB.Select(g => g.Split(',')[0]).Distinct().ToList();
            foreach (string release in releases)
            {
                // read txt file into obsv collection
                var species = genomeDB.Select(g => g.Split(',')[1]).Distinct().ToList();
                Dictionary<string, string> genomes = new Dictionary<string, string>();
                Dictionary<string, string> organisms = new Dictionary<string, string>();

                foreach (string genome in genomeDB.Where(g => g.Contains(release)))
                {
                    var splt = genome.Split(',');
                    genomes.Add(splt[1], splt[3]); // <Species, GenomeVer>
                    organisms.Add(splt[1], splt[2]); // <Species, OrganismName>
                }

                if (genomes.Count > 0)
                {
                    ensemblReleases.Add(new EnsemblRelease() { Release = release, Species = new ObservableCollection<string>(species), Genomes = genomes, Organisms = organisms });
                }
            }
            return ensemblReleases;
        }
    }
}