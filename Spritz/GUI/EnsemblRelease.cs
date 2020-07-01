using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.IO;

namespace SpritzGUI
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
            string releasefolder = Path.Combine(Directory.GetCurrentDirectory(), "EnsemblReleases");
            var releases = Directory.GetFiles(releasefolder, "*.txt").Select(Path.GetFileNameWithoutExtension).ToList();

            var genomeDB = File.ReadAllLines(Path.Combine(Directory.GetCurrentDirectory(), "genomes.csv")).Where(line => !line.StartsWith("#")).ToList();
            foreach (string release in releases)
            {
                // read txt file into obsv collection
                var species = File.ReadAllLines(Path.Combine(releasefolder, $"{release}.txt")).Where(line => !line.StartsWith("#")).ToList();
                Dictionary<string, string> genomes = new Dictionary<string, string>();
                Dictionary<string, string> organisms = new Dictionary<string, string>();

                var unsupported = new List<string>();
                foreach (string genome in genomeDB.Where(g => g.Contains(release)))
                {
                    var splt = genome.Split(',');
                    genomes.Add(splt[1], splt[3]); // <Species, GenomeVer>
                    organisms.Add(splt[1], splt[2]); // <Species, OrganismName>

                    if (!string.Equals(splt[4], "86")) // only add species supported in snpeff (ver 86 ensembl)
                    {
                        unsupported.Add(splt[1]);
                    }
                }

                var supported = species.Where(s => !unsupported.Contains(s)).ToList();
                if (genomes.Count > 0)
                {
                    ensemblReleases.Add(new EnsemblRelease() { Release = release, Species = new ObservableCollection<string>(supported), Genomes = genomes, Organisms = organisms });
                }
            }
            return ensemblReleases;
        }
    }
}