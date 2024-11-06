import { useState, useEffect } from "react";
import "../App.css";
import { CopyBlock, dracula } from "react-code-blocks";

interface Job {
  JobID: string;
  Status: string;
  OutputFiles: string[];
  Timestamp: string;
}

const TestQueue: React.FC = () => {
  const [jobs, setJobs] = useState<Job[]>([]);
  const [error, setError] = useState<string | null>(null);
  const [selectedJob, setSelectedJob] = useState<Job | null>(null);

  const codeText = selectedJob
    ? `import zarrtraj\nimport MDAnalysis as mda\nu = mda.Universe("/path/to/pdb",\n "${selectedJob.OutputFiles[0]}")`
    : "";

  useEffect(() => {
    const fetchJobs = async () => {
      try {
        const response = await fetch(
          "https://rkszuevwgb.execute-api.us-east-1.amazonaws.com/query"
        ); // Adjust the URL to match your API endpoint
        if (!response.ok) {
          throw new Error(`Failed to fetch jobs: ${response.statusText}`);
        }
        const data = await response.json();
        setJobs(data); // Assuming data is an array of jobs
      } catch (err) {
        setError("Failed to fetch jobs. Please try again later.");
        console.error(err);
      }
    };
    fetchJobs();
  }, []);

  const handleInfoClick = (job: Job) => {
    setSelectedJob(job);
  };

  // Function to close the popup
  const closePopup = () => {
    setSelectedJob(null);
  };

  return (
    <div className="submit-jobs-container">
      <h2>Job List</h2>
      {error && <p className="error-message">{error}</p>}
      <div className="scrollable-panel">
        {jobs.map((job) => (
          <div key={job.JobID} className="job-tile">
            <h3>{job.JobID}</h3>
            <p>Status: {job.Status}</p>
            <p>Created: {new Date(job.Timestamp).toLocaleString()}</p>
            <ul>
              {job.OutputFiles.map((file, index) => (
                <li key={index}>
                  <a href={file} target="_blank" rel="noopener noreferrer">
                    {file}
                  </a>
                </li>
              ))}
            </ul>
            <button
              className="see-info-button"
              onClick={() => handleInfoClick(job)}
            >
              info
            </button>
          </div>
        ))}
      </div>
      {selectedJob && (
        <div className="popup-overlay" onClick={closePopup}>
          <div className="popup-content" onClick={(e) => e.stopPropagation()}>
            <h3>Job Info for {selectedJob.JobID}</h3>
            <p>Output files:</p>
            <ul>
              {selectedJob.OutputFiles.map((file, index) => (
                <li key={index}>
                  <a href={file} target="_blank" rel="noopener noreferrer">
                    {file}
                  </a>
                </li>
              ))}
            </ul>
            <p>To access the resulting trajectories, do:</p>
            <div className="code-block">
              <CopyBlock
                text={codeText}
                language="python"
                showLineNumbers={true}
                theme={dracula}
              />
            </div>

            <button className="close-button" onClick={closePopup}>
              ok
            </button>
          </div>
        </div>
      )}
    </div>
  );
};

export default TestQueue;
